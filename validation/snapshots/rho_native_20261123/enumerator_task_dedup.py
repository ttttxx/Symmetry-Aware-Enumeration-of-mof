"""Experimental fixed-point overlap pruning; production entry points unchanged.

Only fully exhausted tasks contribute exclusion regions. A previous partition
P describes equal-color constraints within each block. Compiling those onto a
current task Q lets DFS discard a prefix when every completion is in Fix(P).
The existing canonical-prefix pruning and leaf Trie are retained independently.
"""

import time

from enumerator import MOFEnumerator


class TaskDedupMOFEnumerator(MOFEnumerator):
    """Preserve task order, excluding previously exhausted fixed-point spaces.

    This experiment targets MOFEnumerator's ordinary schedule, not the exact
    rho-quota subclass: a task filtered to one stabilizer order does not cover
    its entire fixed-point domain.
    """

    def _reset_task_schedule(self):
        self._covered_partitions = []
        self.task_dedup_stats = {
            "tasks": 0, "completed_tasks": 0, "pruned_prefixes": 0,
            "leaf_visits": 0, "duplicate_leaves": 0, "compile_s": 0.0,
        }

    def enumerate_subgroups(self, atom_types):
        self._reset_task_schedule()
        return super().enumerate_subgroups(atom_types)

    def _choices(self, orbit, atom_types):
        if self.node_degrees:
            degree = self.node_degrees[orbit[0]]
            if self.topology_requirements and self.user_selected_ligands:
                return self.user_selected_ligands.get(degree, [])
            if degree in self.user_selected_ligands:
                return self.user_selected_ligands[degree]
        return atom_types

    def _compile_exclusions(self, orbits, atom_types):
        """Build equality events and masks of regions decided at each depth.

        An exclusion stays active until an equality is violated. Once all its
        equalities are decided, an active exclusion covers every completion.
        Singleton-label domains are treated as constants even before assignment.
        """
        size = len(orbits)
        node_orbit = [0] * self.n
        for index, orbit in enumerate(orbits):
            for node in orbit:
                node_orbit[node] = index
        domains = [tuple(self._choices(orbit, atom_types)) for orbit in orbits]
        if any(not domain for domain in domains):
            return [[] for _ in orbits], [1] * (size + 1), 1
        grouped_events = [{} for _ in orbits]
        decided = [0] * (size + 1)
        all_active = 0
        for partition in self._covered_partitions:
            parent = list(range(size))

            def find(index):
                while parent[index] != index:
                    parent[index] = parent[parent[index]]
                    index = parent[index]
                return index

            for block in partition:
                anchor = node_orbit[block[0]]
                for node in block[1:]:
                    parent[find(node_orbit[node])] = find(anchor)
            components = {}
            for index in range(size):
                components.setdefault(find(index), []).append(index)
            constraints = []
            impossible = False
            for component in components.values():
                if len(component) < 2:
                    continue
                constants = {domains[i][0] for i in component if len(domains[i]) == 1}
                if len(constants) > 1:
                    impossible = True
                    break
                variables = [i for i in component if len(domains[i]) > 1]
                if constants:
                    constant = next(iter(constants))
                    if any(constant not in domains[i] for i in variables):
                        impossible = True
                        break
                    constraints.extend((i, "constant", constant) for i in variables)
                elif variables:
                    if not set.intersection(*(set(domains[i]) for i in variables)):
                        impossible = True
                        break
                    anchor = variables[0]
                    constraints.extend((i, "orbit", anchor) for i in variables[1:])
            if impossible:
                continue  # This earlier region cannot intersect this task.
            bit = all_active + 1  # Active bits are allocated consecutively.
            all_active |= bit
            last_depth = 0
            for index, kind, other in constraints:
                event = (kind, other)
                grouped_events[index][event] = grouped_events[index].get(event, 0) | bit
                last_depth = max(last_depth, index + 1)
            decided[last_depth] |= bit
        for depth in range(1, size + 1):
            decided[depth] |= decided[depth - 1]
        events = [list(group.items()) for group in grouped_events]
        return events, decided, all_active

    def _try_add_result(self, coloring):
        self.task_dedup_stats["leaf_visits"] += 1
        before = len(self.classes)
        result = super()._try_add_result(coloring)
        if len(self.classes) == before:
            self.task_dedup_stats["duplicate_leaves"] += 1
        return result

    def _search_compiled_task(self, orbits, atom_types, current_coloring,
                              active_perm_indices, current_degree_counts,
                              identity_start, identity_budget, deadline):
        """Hook for a fused implementation without per-level method wrapping."""
        return super()._dfs_lattice_pruning(
            0, orbits, atom_types, current_coloring, active_perm_indices,
            current_degree_counts, identity_start, identity_budget, deadline)

    def _dfs_lattice_pruning(self, orbit_idx, orbits, atom_types,
                             current_coloring, active_perm_indices,
                             current_degree_counts=None, identity_start=None,
                             identity_budget=None, deadline=None):
        # Keep direct calls to the inherited DFS usable without a task schedule.
        if not hasattr(self, "_covered_partitions"):
            return super()._dfs_lattice_pruning(
                orbit_idx, orbits, atom_types, current_coloring, active_perm_indices,
                current_degree_counts, identity_start, identity_budget, deadline)
        if orbit_idx == 0:
            self.task_dedup_stats["tasks"] += 1
            started = time.perf_counter()
            self._exclusion_events, self._decided_masks, initial = self._compile_exclusions(orbits, atom_types)
            self.task_dedup_stats["compile_s"] += time.perf_counter() - started
            self._exclusion_stack = [initial] + [0] * len(orbits)
            mask = initial
        else:
            mask = self._exclusion_stack[orbit_idx - 1]
            if mask:
                current_value = current_coloring[orbits[orbit_idx - 1][0]]
                for (kind, other), exclusions in self._exclusion_events[orbit_idx - 1]:
                    if not (mask & exclusions):
                        continue
                    other_value = other if kind == "constant" else current_coloring[orbits[other][0]]
                    if current_value != other_value:
                        mask &= ~exclusions
                        if not mask:
                            break
            self._exclusion_stack[orbit_idx] = mask
        if mask & self._decided_masks[orbit_idx]:
            self.task_dedup_stats["pruned_prefixes"] += 1
            result = None
        elif orbit_idx == 0:
            result = self._search_compiled_task(
                orbits, atom_types, current_coloring, active_perm_indices,
                current_degree_counts, identity_start, identity_budget, deadline)
        else:
            result = super()._dfs_lattice_pruning(
                orbit_idx, orbits, atom_types, current_coloring, active_perm_indices,
                current_degree_counts, identity_start, identity_budget, deadline)
        if orbit_idx == 0:
            truncated = (
                (deadline is not None and time.monotonic() >= deadline)
                or (self.max_total_results is not None and len(self.classes) >= self.max_total_results)
                or (identity_start is not None and identity_budget is not None
                    and len(self.classes) - identity_start >= identity_budget)
            )
            if not truncated:
                self._covered_partitions.append(tuple(tuple(orbit) for orbit in orbits))
                self.task_dedup_stats["completed_tasks"] += 1
        return result
