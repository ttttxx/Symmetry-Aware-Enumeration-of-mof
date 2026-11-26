"""Opt-in exact-stabilizer-layer quota experiment.

The production implementation remains in :mod:`enumerator`.  This subclass
reuses its canonicalization and DFS pruning, but schedules fixed-point tasks by
subgroup order and retains at most a configured number of representatives with
each exact stabilizer order.  Within a layer, the inherited deterministic DFS
order is used; this is a reproducible bounded prefix, not a statistical sample.
"""

import time
from typing import Dict, List, Tuple

from enumerator import MOFEnumerator, _TrieNode


class RhoQuotaMOFEnumerator(MOFEnumerator):
    def __init__(self, *args, rho_quota_per_level: int, **kwargs):
        if rho_quota_per_level <= 0:
            raise ValueError("rho_quota_per_level must be a positive integer")
        super().__init__(*args, **kwargs)
        self.rho_quota_per_level = int(rho_quota_per_level)
        self.rho_layer_report: List[Dict] = []
        self.rho_setup_time_s = 0.0
        self._active_rho_order = None
        self._active_rho_count = 0
        self._rho_seen = set()

    def _try_add_result(self, coloring: List[str]):
        canonical = self._get_canonical(tuple(coloring))
        if canonical in self._rho_seen or self._trie_lookup(canonical) is not None:
            return
        self._rho_seen.add(canonical)

        stabilizer_size, _ = self._compute_stabilizer_and_orbit(canonical)
        if stabilizer_size != self._active_rho_order:
            return

        class_id = len(self.classes)
        self._trie_insert(canonical, class_id, stabilizer_size)
        self.classes[class_id] = {
            "id": class_id,
            "coloring": list(canonical),
            "stabilizer_size": stabilizer_size,
            "orbit_size": len(self.perms) // stabilizer_size,
        }
        self._active_rho_count += 1

    def _dfs_lattice_pruning(self, *args, **kwargs):
        if self._active_rho_count >= self.rho_quota_per_level:
            return
        if self.max_total_results is not None and len(self.classes) >= self.max_total_results:
            return
        return super()._dfs_lattice_pruning(*args, **kwargs)

    def enumerate_subgroups(self, atom_types: List[str]):
        """Enumerate exact stabilizer-order layers, highest order first."""
        try:
            from sympy.combinatorics import Permutation, PermutationGroup
        except Exception as exc:
            raise RuntimeError(
                "The rho-quota experiment requires SymPy's complete subgroup enumeration"
            ) from exc

        self._last_atom_types = atom_types
        self.classes = {}
        self._trie_root = _TrieNode()
        self._rho_seen = set()
        self.rho_layer_report = []
        setup_started_at = time.monotonic()

        def subgroup_partitions(group):
            """Return all subgroups using the same closure construction as the legacy enumerator."""
            group_order = group.order()
            found = []
            seen = set()

            def add_if_new(subgroup):
                key = frozenset(tuple(elem.array_form) for elem in subgroup.generate())
                if key in seen:
                    return False
                seen.add(key)
                found.append(subgroup)
                return True

            for element in group.generate():
                add_if_new(PermutationGroup(element))

            index = 0
            while index < len(found):
                current = found[index]
                for earlier in found[:index]:
                    combined = PermutationGroup(
                        *(list(current.generators) + list(earlier.generators))
                    )
                    if combined.order() < group_order:
                        add_if_new(combined)
                index += 1
            return found

        def to_tuple(permutation):
            values = list(permutation.array_form)
            values += list(range(len(values), self.n))
            return tuple(values)

        group_permutations = [
            Permutation(list(permutation), size=self.n) for permutation in self.perms
        ]
        group = PermutationGroup(*group_permutations)
        group_order = group.order()
        if group_order != len(self.perms):
            raise RuntimeError(
                f"Permutation list has {len(self.perms)} entries, but generated group "
                f"has order {group_order}; exact rho layers are undefined for this input"
            )

        # Keep partitions grouped by subgroup order.  A partition may occur at
        # several orders; retaining each order is safe because leaves are
        # filtered by their actual full-group stabilizer.
        layers: Dict[int, Dict[Tuple[Tuple[int, ...], ...], List[List[int]]]] = {}

        def add_partition(order, orbits):
            fingerprint = tuple(tuple(orbit) for orbit in orbits)
            layers.setdefault(order, {})[fingerprint] = orbits

        add_partition(group_order, self._compute_orbits_by_gen(self.perms))
        for subgroup in subgroup_partitions(group):
            order = subgroup.order()
            if order <= 1 or order >= group_order:
                continue
            generators = [to_tuple(generator) for generator in subgroup.generators]
            add_partition(order, self._compute_orbits_by_gen(generators))

        identity_orbits = [[index] for index in range(self.n)]
        add_partition(1, identity_orbits)
        self.rho_setup_time_s = time.monotonic() - setup_started_at

        ordered_layers = sorted(layers, reverse=True)
        active_permutations = list(range(1, len(self.perms)))
        stopped_by_global_cap = False

        for layer_position, stabilizer_order in enumerate(ordered_layers):
            layer_started_at = time.monotonic()
            layer_tasks = sorted(layers[stabilizer_order].items(), key=lambda item: item[0])
            self._active_rho_order = stabilizer_order
            self._active_rho_count = 0
            tasks_started = 0
            task_timeouts = 0

            for _, orbits in layer_tasks:
                if self._active_rho_count >= self.rho_quota_per_level:
                    break
                if self.max_total_results is not None and len(self.classes) >= self.max_total_results:
                    stopped_by_global_cap = True
                    break

                tasks_started += 1
                deadline = (
                    time.monotonic() + self.task_time_limit
                    if self.task_time_limit else None
                )
                self._dfs_lattice_pruning(
                    0,
                    orbits,
                    atom_types,
                    [None] * self.n,
                    active_permutations,
                    None,
                    None,
                    deadline,
                )
                if (
                    self.max_total_results is not None
                    and len(self.classes) >= self.max_total_results
                ):
                    # The DFS has no completion flag. Treat reaching the cap
                    # as truncation even if the final possible leaf happened
                    # to coincide with the cap; this avoids claiming that a
                    # fixed layer was exhausted without evidence.
                    stopped_by_global_cap = True
                    break
                if deadline is not None and time.monotonic() >= deadline:
                    task_timeouts += 1

            if stopped_by_global_cap:
                status = "global_output_cap"
            elif self._active_rho_count >= self.rho_quota_per_level:
                status = "quota_reached"
            elif task_timeouts:
                status = "time_limited"
            else:
                status = "layer_exhausted"

            report = {
                "stabilizer_order": stabilizer_order,
                "rho": stabilizer_order / group_order,
                "quota": self.rho_quota_per_level,
                "returned": self._active_rho_count,
                "fixed_layer_tasks": len(layer_tasks),
                "tasks_started": tasks_started,
                "task_timeouts": task_timeouts,
                "elapsed_s": time.monotonic() - layer_started_at,
                "status": status,
            }
            self.rho_layer_report.append(report)
            if self.verbose:
                print(
                    f"    [rho layer] |Stab|={stabilizer_order} "
                    f"(rho={stabilizer_order / group_order:.6g}): "
                    f"returned={self._active_rho_count}/{self.rho_quota_per_level}, "
                    f"tasks={tasks_started}/{len(layer_tasks)}, status={status}"
                )

            if stopped_by_global_cap:
                for unvisited_order in ordered_layers[layer_position + 1:]:
                    self.rho_layer_report.append({
                        "stabilizer_order": unvisited_order,
                        "rho": unvisited_order / group_order,
                        "quota": self.rho_quota_per_level,
                        "returned": 0,
                        "fixed_layer_tasks": len(layers[unvisited_order]),
                        "tasks_started": 0,
                        "task_timeouts": 0,
                        "status": "not_started_global_output_cap",
                    })
                break

        self._active_rho_order = None
