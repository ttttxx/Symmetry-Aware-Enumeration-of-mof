"""Exhaustive comparison experiments: product counter vs canonical-prefix DFS.

itertools.product generates full orbit-label tuples in C; it does not perform
DFS prefix pruning. All variants retain integer encoding, the same group,
Trie storage and full representative metadata. Single-task variants omit the
subgroup priority schedule, so bounded outputs can differ from production.
"""

import itertools
import math
import time

from enumerator_casm_inspired import CasmInspiredMOFEnumerator


class CounterMOFEnumerator(CasmInspiredMOFEnumerator):
    """Keep subgroup tasks, replace each DFS by full-candidate iteration."""

    def _search_compiled_task(self, orbits, atom_types, coloring,
                              active_perm_indices, current_degree_counts,
                              identity_start, identity_budget, deadline):
        n = self.n
        domains = [self._choices(orbit, atom_types) for orbit in orbits]
        identity = tuple(range(n))
        inverses = [inverse for inverse in self.perms_inv if inverse != identity]
        if self.perms.count(identity) != 1:
            raise ValueError("Counter experiment requires exactly one identity operation")
        node_to_orbit = [0] * n
        for index, orbit in enumerate(orbits):
            for node in orbit:
                node_to_orbit[node] = index
        direct = len(orbits) == n and all(orbit == [i] for i, orbit in enumerate(orbits))
        # Compare views of the orbit-label tuple. Remove positions guaranteed
        # equal by this task or by fixed single-label domains, while preserving
        # the original site order of all potentially differing positions.
        comparison_pairs = []
        for inverse in inverses:
            pairs = []
            for site in range(n):
                original_var, image_var = node_to_orbit[site], node_to_orbit[inverse[site]]
                if original_var == image_var:
                    continue
                a, b = domains[original_var], domains[image_var]
                if len(a) == len(b) == 1 and a[0] == b[0]:
                    continue
                pairs.append((original_var, image_var))
            comparison_pairs.append(pairs)
        events = [(i, event) for i, event in enumerate(self._exclusion_events) if event]
        decided = self._decided_masks
        initial_mask = self._exclusion_stack[0]
        stats = self.task_dedup_stats
        stats.setdefault("counter_candidates", 0)
        stats.setdefault("counter_excluded", 0)
        # All points of a task assign the same number of sites of each degree.
        if self.topology_requirements:
            counts = dict(current_degree_counts or {degree: 0 for degree in self.topology_requirements})
            for orbit in orbits:
                if self.node_degrees:
                    degree = self.node_degrees[orbit[0]]
                    counts[degree] += len(orbit)
                    if counts[degree] > self.topology_requirements.get(degree, counts[degree]):
                        return
        for assignment in itertools.product(*domains):
            if deadline is not None and time.monotonic() >= deadline:
                break
            if self.max_total_results is not None and len(self.classes) >= self.max_total_results:
                break
            if identity_start is not None and identity_budget is not None:
                if len(self.classes) - identity_start >= identity_budget:
                    break
            stats["counter_candidates"] += 1
            mask = initial_mask
            excluded = False
            if mask:
                for depth, comparisons in events:
                    value = assignment[depth]
                    for (kind, other), regions in comparisons:
                        if mask & regions:
                            other_value = other if kind == "constant" else assignment[other]
                            if value != other_value:
                                mask &= ~regions
                    if mask & decided[depth + 1]:
                        excluded = True
                        break
                    if not mask:
                        break
            if excluded:
                stats["counter_excluded"] += 1
                continue
            stabilizer = 1
            canonical = True
            for pairs in comparison_pairs:
                for i, j in pairs:
                    original, image = assignment[i], assignment[j]
                    if image < original:
                        canonical = False
                        break
                    if image > original:
                        break
                else:
                    stabilizer += 1
                if not canonical:
                    break
            if canonical:
                candidate = assignment if direct else tuple(assignment[i] for i in node_to_orbit)
                self._store_canonical_leaf(candidate, stabilizer)


class _SingleTaskSchedule:
    """One identity task, to isolate traversal from subgroup scheduling."""

    def enumerate_subgroups(self, atom_types):
        self._reset_task_schedule()
        self._last_atom_types = atom_types
        orbits = [[i] for i in range(self.n)]
        naive = math.prod(len(self._choices(orbit, atom_types)) for orbit in orbits)
        budget = None
        if self._identity_orbit_limit is not None and naive / max(1, len(self.perms)) > self._identity_orbit_limit:
            budget = self._identity_new_limit
        deadline = time.monotonic() + self.task_time_limit if self.task_time_limit else None
        self._dfs_lattice_pruning(
            0, orbits, atom_types, [None] * self.n, list(range(1, len(self.perms))),
            identity_start=0 if budget is not None else None,
            identity_budget=budget, deadline=deadline)


class SingleTaskDfsMOFEnumerator(_SingleTaskSchedule, CasmInspiredMOFEnumerator):
    """Optimized DFS without subgroup preparation or priority tasks."""


class SingleTaskCounterMOFEnumerator(_SingleTaskSchedule, CounterMOFEnumerator):
    """CASM-like full-space counter plus lazy canonicality filtering."""
