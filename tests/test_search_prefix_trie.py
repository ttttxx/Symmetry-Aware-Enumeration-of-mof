"""Completed-prefix certificates must save search without losing completions."""
import contextlib
import io
import itertools

import pytest

from enumerator import MOFEnumerator
from enumerator_casm_inspired import CasmInspiredMOFEnumerator
from enumerator_storage import HashStorageMOFEnumerator
from enumerator_native import HybridMOFEnumerator, native_available
from search_prefix_trie import CompletedPrefixCache
from tests.test_subgroup_task_dedup import dihedral, run


def interrupted_schedule(self, labels):
    """First fixed-point task ends at a deterministic new-result budget."""
    source = [[0, 1]] + [[i] for i in range(2, self.n)]
    target = [[i] for i in range(self.n)]
    for partition, budget in ((source, self.test_budget), (target, None)):
        self._dfs_lattice_pruning(0, partition, labels, [None] * self.n,
                                 list(range(1, len(self.perms))),
                                 identity_start=len(self.classes) if budget is not None else None,
                                 identity_budget=budget)


def run_interrupted(backend, enabled, budget=4, n=4, max_nodes=50000):
    identity = tuple(range(n))
    swap = (1, 0) + tuple(range(2, n))
    instance = HybridMOFEnumerator([identity, swap], n, backend=backend,
                                   prefix_trie=enabled, prefix_trie_max_nodes=max_nodes,
                                   verbose=False, max_total_results=None,
                                   identity_orbit_limit=float("inf"))
    instance.test_budget = budget
    with contextlib.redirect_stdout(io.StringIO()):
        records = instance.enumerate(["A", "B"])
    return instance, records


def test_cache_requires_closed_prefix_equalities_and_domain_containment():
    cache = CompletedPrefixCache()
    source = [[0, 1], [2], [3]]
    domains = [[0, 1]] * 3
    # Entire source prefix 0 finished; source prefix 1 / 0 is still open.
    cache.capture(source, domains, [1, 0], [1, 0, 0], 2)
    query = cache.task_query([[0], [1], [2], [3]], [[0, 1]] * 4)
    assert not query.covers([0, None, None, None], 1)  # Equality not forced.
    assert query.covers([0, 0, None, None], 2)       # Skips four leaves.
    assert not query.covers([1, 1, 0, None], 3)     # Existing open path.
    assert not query.covers([0, 1, 0, 0], 4)
    wider = cache.task_query([[i] for i in range(4)], [[0, 1, 2]] * 4)
    assert not wider.covers([0, 0, None, None], 2)  # Source lacked label 2.
    # Reordered variables and singleton domains may prove equality at root.
    reordered = cache.task_query([[3], [1], [0], [2]], [[0, 1], [0], [0], [0, 1]])
    assert reordered.covers([None] * 4, 0)
    assert cache.stats["early_hits"] == 2


@pytest.mark.parametrize("backend", ["python", "cpp"])
@pytest.mark.parametrize("budget", range(9))
@pytest.mark.parametrize("max_nodes", [0, 2, 50000])
def test_every_interruption_frontier_preserves_complete_records(monkeypatch, backend, budget, max_nodes):
    if backend == "cpp" and not native_available():
        pytest.skip("build optional C++ core")
    monkeypatch.setattr(MOFEnumerator, "enumerate_subgroups", interrupted_schedule)
    baseline, expected = run_interrupted(backend, False, budget)
    instance, records = run_interrupted(backend, True, budget, max_nodes=max_nodes)
    assert records == expected
    # Independent orbit/metadata oracle for the two-element swap group.
    canonical = {min(c, (c[1], c[0], c[2], c[3]))
                 for c in itertools.product("AB", repeat=4)}
    assert {tuple(r["coloring"]) for r in records} == canonical
    assert sum(r["orbit_size"] for r in records) == 16
    assert instance.prefix_trie_stats["peak_nodes"] <= max_nodes
    assert instance.task_dedup_stats["search_nodes"] <= baseline.task_dedup_stats["search_nodes"]
    if budget == 4 and max_nodes == 50000:
        assert instance.prefix_trie_stats["early_hits"] == 1
        assert instance.task_dedup_stats["duplicate_leaves"] == 0
        assert baseline.task_dedup_stats["duplicate_leaves"] == 4
        assert instance.task_dedup_stats["search_nodes"] < baseline.task_dedup_stats["search_nodes"]


@pytest.mark.parametrize("cls,options", [(CasmInspiredMOFEnumerator, {}),
                                         (HashStorageMOFEnumerator, {}),
                                         (HybridMOFEnumerator, {"backend": "cpp"})])
@pytest.mark.parametrize("limits", [{}, {"max_total_results": 7},
                                    {"identity_orbit_limit": 0, "identity_new_limit": 2}])
def test_existing_storage_schedule_and_budgets_unchanged(cls, options, limits):
    if options.get("backend") == "cpp" and not native_available():
        pytest.skip("build optional C++ core")
    _, expected = run(cls, dihedral(5), **options, **limits)
    instance, records = run(cls, dihedral(5), prefix_trie=True, **options, **limits)
    assert records == expected
    if not limits:
        assert instance.prefix_trie_stats["nodes_created"] == 0
    with contextlib.redirect_stdout(io.StringIO()):
        repeated = instance.enumerate(["A", "B", "C"])
    assert repeated == expected


def test_counter_rejects_prefix_mode():
    with pytest.raises(ValueError, match="DFS|dfs"):
        HybridMOFEnumerator([(0,)], 1, prefix_trie=True, native_kernel="counter")


def test_cli_prefix_mode_is_explicit_and_reports_configuration(monkeypatch, capsys):
    import main
    cls, options = main.get_enumerator_implementation("cpp", True, 123)
    assert cls is HybridMOFEnumerator
    assert options == dict(backend="cpp", prefix_trie=True, prefix_trie_max_nodes=123)
    with pytest.raises(ValueError, match="requires"):
        main.get_enumerator_implementation("legacy", True)
    monkeypatch.setattr("sys.argv", ["main.py", "--enumeration-backend", "legacy", "--prefix-trie"])
    with pytest.raises(SystemExit) as exc:
        main.main()
    assert exc.value.code == 2
    assert "requires --enumeration-backend" in capsys.readouterr().err
    class InputReached(Exception):
        pass
    def stop_before_input(*args, **kwargs):
        raise InputReached
    monkeypatch.setattr(main, "load_graph", stop_before_input)
    monkeypatch.setattr("sys.argv", ["main.py", "--prefix-trie"])
    with pytest.raises(InputReached):
        main.main()


@pytest.mark.parametrize("backend", ["python", "cpp"])
def test_expired_prefix_task_never_claims_a_subtree(backend):
    if backend == "cpp" and not native_available():
        pytest.skip("build optional C++ core")
    # Native normal schedule also exercises a deadline that expires on entry.
    instance = HybridMOFEnumerator(dihedral(4), 4, backend=backend, prefix_trie=True,
                                   verbose=False, max_total_results=None,
                                   task_time_limit=1e-12)
    with contextlib.redirect_stdout(io.StringIO()):
        result = instance.enumerate(["A", "B"])
    assert not result
    assert not instance._covered_partitions
    assert instance.prefix_trie_stats["closed_prefixes"] == 0


@pytest.mark.parametrize("backend", ["python", "cpp"])
@pytest.mark.parametrize("budget", [1, 3, 5, 7])
def test_multiple_interrupted_non_nested_tasks_do_not_lose_canonical_leaves(monkeypatch, backend, budget):
    if backend == "cpp" and not native_available():
        pytest.skip("build optional C++ core")
    def schedule(self, labels):
        for partition in ([[[0, 2], [1, 3]], [[0, 1], [2, 3]], [[0, 3], [1, 2]]]):
            self._dfs_lattice_pruning(0, partition, labels, [None] * 4,
                                     list(range(1, len(self.perms))),
                                     identity_start=len(self.classes), identity_budget=budget)
        self._dfs_lattice_pruning(0, [[i] for i in range(4)], labels, [None] * 4,
                                 list(range(1, len(self.perms))))
    monkeypatch.setattr(MOFEnumerator, "enumerate_subgroups", schedule)
    _, expected = run(HybridMOFEnumerator, dihedral(4), backend=backend)
    instance, records = run(HybridMOFEnumerator, dihedral(4), backend=backend, prefix_trie=True)
    assert records == expected
    canonical = {min(tuple(c[p.index(i)] for i in range(4)) for p in dihedral(4))
                 for c in itertools.product("ABC", repeat=4)}
    assert {tuple(r["coloring"]) for r in records} == canonical
    assert sum(r["orbit_size"] for r in records) == 81
    if budget <= 5:
        assert instance.prefix_trie_stats["snapshots"] >= 1


def test_partial_prefix_reuse_also_retains_full_representative_trie(monkeypatch):
    monkeypatch.setattr(MOFEnumerator, "enumerate_subgroups", interrupted_schedule)
    options = dict(verbose=False, max_total_results=None, identity_orbit_limit=float("inf"))
    def execute(enabled):
        instance = CasmInspiredMOFEnumerator([(0, 1, 2, 3), (1, 0, 2, 3)], 4,
                                             prefix_trie=enabled, **options)
        instance.test_budget = 4
        with contextlib.redirect_stdout(io.StringIO()):
            records = instance.enumerate(["A", "B"])
        return instance, records
    _, expected = execute(False)
    instance, records = execute(True)
    assert records == expected
    assert instance.prefix_trie_stats["early_hits"] == 1
    assert not hasattr(instance, "_canonical_index")
    assert all(instance._trie_lookup(tuple(r["coloring"])) is not None for r in records)
    # Reuse with different label codes must rebuild the prefix certificates.
    with contextlib.redirect_stdout(io.StringIO()):
        changed = instance.enumerate(["Y", "Z"])
    assert len(changed) == 12
    assert instance.prefix_trie_stats["snapshots"] == 1
    assert instance.prefix_trie_stats["early_hits"] == 1
