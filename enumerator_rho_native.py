"""Exact-stabilizer quotas with integer DFS, native search and layer-scoped coverage."""

from enumerator_native import HybridMOFEnumerator
from enumerator_rho_quota import RhoQuotaScheduleMixin


class HybridRhoQuotaMOFEnumerator(RhoQuotaScheduleMixin, HybridMOFEnumerator):
    """Keep the quota task/label order, admitting only the current exact order.

    Complete-task and completed-prefix certificates describe target-order search,
    so their lifetime is one stabilizer layer. Output keys remain global.
    """

    def _run_schedule(self, atom_types):
        self._canonical_index = {}
        self._use_bytes_key = self.storage_key == "bytes" and len(self._code_to_label) <= 256
        self._reset_task_schedule()
        self._coverage_order = None
        try:
            return self._enumerate_quota_subgroups(atom_types)
        finally:
            self._active_rho_order = None

    def _store_canonical_leaf(self, canonical, stabilizer_size):
        if stabilizer_size != self._active_rho_order:
            self.task_dedup_stats["leaf_visits"] += 1
            return
        return super()._store_canonical_leaf(canonical, stabilizer_size)

    def _run_rho_task(self, orbits, atom_types, active_permutations, deadline):
        if self._coverage_order != self._active_rho_order:
            self._covered_partitions = []
            if self.prefix_trie:
                from search_prefix_trie import CompletedPrefixCache
                self._prefix_cache = CompletedPrefixCache(self.prefix_trie_max_nodes)
                self.prefix_trie_stats = self._prefix_cache.stats
            self._coverage_order = self._active_rho_order
            # Native certificates reset at the same boundary, before any task
            # can be skipped by the Python coverage check.
            if self._native_context is not None:
                if self._core.mof_core_target_order(self._native_context, self._active_rho_order):
                    raise RuntimeError(self._core.mof_core_error().decode())
        before = len(self.classes)
        self._dfs_lattice_pruning(
            0, orbits, atom_types, [None] * self.n, active_permutations,
            identity_start=before,
            identity_budget=self.rho_quota_per_level - self._active_rho_count,
            deadline=deadline,
        )
        self._active_rho_count += len(self.classes) - before
