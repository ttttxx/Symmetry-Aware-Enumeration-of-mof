"""Experimental integer DFS with overlap exclusion and lazy group comparison.

Retains the ordinary subgroup schedule, output ordering and leaf Trie. Labels
are integer codes internally and decoded before returning. Each permutation's
known-equal prefix is carried down DFS; surviving operations at a complete
canonical leaf give its stabilizer directly, without materializing its orbit.
"""

import time

from enumerator_task_dedup import TaskDedupMOFEnumerator


class CasmInspiredMOFEnumerator(TaskDedupMOFEnumerator):
    def enumerate(self, atom_types, composition=None):
        original_pools = self.user_selected_ligands
        labels = sorted(set(atom_types).union(*(set(pool) for pool in original_pools.values())))
        self._code_to_label = labels
        self._label_to_code = {label: i for i, label in enumerate(labels)}
        self._labels_encoded = True
        self.user_selected_ligands = {
            degree: [self._label_to_code[label] for label in pool]
            for degree, pool in original_pools.items()
        }
        try:
            return super().enumerate([self._label_to_code[label] for label in atom_types], composition)
        finally:
            for record in self.classes.values():
                record["coloring"] = [labels[value] for value in record["coloring"]]
            self.user_selected_ligands = original_pools
            self._last_atom_types = atom_types
            self._labels_encoded = False

    def _trie_lookup(self, canonical):
        if hasattr(self, "_label_to_code") and not self._labels_encoded:
            try:
                canonical = tuple(self._label_to_code[label] for label in canonical)
            except KeyError:
                return None
        return super()._trie_lookup(canonical)

    def _get_canonical(self, coloring):
        # Fallback for explicit leaf calls: compare permutation views, allocate
        # only the final winner. Normal fused DFS already certifies its leaves.
        best = tuple(range(self.n))
        for inverse in self.perms_inv:
            for index in range(self.n):
                candidate, incumbent = coloring[inverse[index]], coloring[best[index]]
                if candidate < incumbent:
                    best = inverse
                    break
                if candidate > incumbent:
                    break
        return tuple(coloring[index] for index in best)

    def _stabilizer_size(self, coloring):
        return sum(all(coloring[i] == coloring[inverse[i]] for i in range(self.n))
                   for inverse in self.perms_inv)

    def _store_canonical_leaf(self, canonical, stabilizer_size):
        self.task_dedup_stats["leaf_visits"] += 1
        if self._trie_lookup(canonical) is not None:
            self.task_dedup_stats["duplicate_leaves"] += 1
            return
        cid = len(self.classes)
        self._trie_insert(canonical, cid, stabilizer_size)
        self.classes[cid] = dict(id=cid, coloring=list(canonical),
                                 stabilizer_size=stabilizer_size,
                                 orbit_size=len(self.perms) // stabilizer_size)

    def _try_add_result(self, coloring):
        canonical = self._get_canonical(tuple(coloring))
        self._store_canonical_leaf(canonical, self._stabilizer_size(canonical))

    def _search_compiled_task(self, orbits, atom_types, coloring,
                              active_perm_indices, current_degree_counts,
                              identity_start, identity_budget, deadline):
        inverse_permutations = self.perms_inv
        identity = tuple(range(self.n))
        if self.perms.count(identity) != 1:
            raise ValueError("The experimental DFS requires exactly one identity operation")
        # Use the full group independently of the caller's identity position.
        active = [(i, 0) for i, permutation in enumerate(self.perms) if permutation != identity]
        domains = [self._choices(orbit, atom_types) for orbit in orbits]
        degrees = [self.node_degrees[orbit[0]] if self.node_degrees else None for orbit in orbits]
        counts = current_degree_counts
        if counts is None and self.topology_requirements:
            counts = {degree: 0 for degree in self.topology_requirements}
        events, decided = self._exclusion_events, self._decided_masks
        stats = self.task_dedup_stats
        total_limit = self.max_total_results
        n, depth_limit = self.n, len(orbits)

        def search(depth, comparisons, mask):
            if deadline is not None and time.monotonic() >= deadline:
                return
            if total_limit is not None and len(self.classes) >= total_limit:
                return
            if identity_start is not None and identity_budget is not None:
                if len(self.classes) - identity_start >= identity_budget:
                    return
            if depth and mask:
                value = coloring[orbits[depth - 1][0]]
                for (kind, other), regions in events[depth - 1]:
                    if mask & regions:
                        other_value = other if kind == "constant" else coloring[orbits[other][0]]
                        if value != other_value:
                            mask &= ~regions
                            if not mask:
                                break
            if mask & decided[depth]:
                stats["pruned_prefixes"] += 1
                return
            if depth == depth_limit:
                # Each surviving nonidentity operation agrees at every site.
                self._store_canonical_leaf(tuple(coloring), len(comparisons) + 1)
                return
            nodes, degree = orbits[depth], degrees[depth]
            for value in domains[depth]:
                if deadline is not None and time.monotonic() >= deadline:
                    break
                for node in nodes:
                    coloring[node] = value
                tracked = bool(self.topology_requirements and counts is not None and degree is not None)
                if tracked:
                    new_count = counts[degree] + len(nodes)
                    if new_count > self.topology_requirements.get(degree, new_count):
                        for node in nodes:
                            coloring[node] = None
                        continue
                    counts[degree] = new_count
                next_comparisons = []
                pruned = False
                for permutation_index, cursor in comparisons:
                    inverse = inverse_permutations[permutation_index]
                    while cursor < n:
                        current = coloring[cursor]
                        image = coloring[inverse[cursor]]
                        if current is None or image is None:
                            next_comparisons.append((permutation_index, cursor))
                            break
                        if image < current:
                            pruned = True
                            break
                        if image > current:
                            break
                        cursor += 1
                    else:
                        next_comparisons.append((permutation_index, n))
                    if pruned:
                        break
                if not pruned:
                    search(depth + 1, next_comparisons, mask)
                if tracked:
                    counts[degree] -= len(nodes)
                for node in nodes:
                    coloring[node] = None

        return search(0, active, self._exclusion_stack[0])
