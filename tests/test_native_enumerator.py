"""Native/Python equivalence and fallback contracts for subgroup scheduling."""
import contextlib
import io
import itertools

import pytest

from tests.test_subgroup_task_dedup import dihedral, run
from enumerator_storage import HashStorageMOFEnumerator
import enumerator_native as native
from enumerator_native import HybridMOFEnumerator
from enumerator_task_dedup import TaskDedupMOFEnumerator


@pytest.fixture
def core():
    if not native.native_available():
        pytest.skip("Build optional core: python tools/build_enumeration_core.py")


@pytest.mark.parametrize("store", ["trie", "hash"])
@pytest.mark.parametrize("kernel", ["dfs", "counter"])
@pytest.mark.parametrize("perms", [dihedral(4), dihedral(5), list(itertools.permutations(range(3)))])
def test_native_matches_complete_python_records(core, kernel, perms, store):
    _, expected = run(HashStorageMOFEnumerator, perms)
    instance, result = run(HybridMOFEnumerator, perms, backend="cpp", native_kernel=kernel, native_store=store)
    assert result == expected
    assert instance.active_backend == "cpp"
    assert instance._native_context is None
    if store == "hash":
        assert len(instance._canonical_index) == len(result)
    else:
        assert not instance._canonical_index
        assert instance.result_index_stats["hash_keys"] == 0
        assert instance.result_index_stats["nodes"] <= 2 * len(result) + 1
    assert all(instance._trie_lookup(tuple(r["coloring"])) is not None for r in result)


@pytest.mark.parametrize("kernel", ["dfs", "counter"])
@pytest.mark.parametrize("limits", [dict(max_total_results=7),
                                   dict(identity_orbit_limit=0, identity_new_limit=2)])
def test_native_preserves_count_budgets(core, kernel, limits):
    _, expected = run(HashStorageMOFEnumerator, dihedral(4), **limits)
    _, result = run(HybridMOFEnumerator, dihedral(4), backend="cpp", native_kernel=kernel, **limits)
    assert result == expected


def test_python_and_missing_library_fallback(monkeypatch):
    _, expected = run(HashStorageMOFEnumerator, dihedral(4))
    def missing():
        raise OSError("test missing optional library")
    monkeypatch.setattr(native, "_load_core", missing)
    instance, result = run(HybridMOFEnumerator, dihedral(4), backend="auto")
    assert instance.active_backend == "python"
    assert result == expected
    with pytest.raises(RuntimeError, match="build"):
        run(HybridMOFEnumerator, dihedral(4), backend="cpp")
    _, forced = run(HybridMOFEnumerator, dihedral(4), backend="python")
    assert forced == expected


def test_large_label_alphabet_falls_back_without_truncation():
    labels = [f"label_{i:03}" for i in range(257)]
    instance = HybridMOFEnumerator([(0,)], 1, backend="auto", verbose=False, max_total_results=None)
    with contextlib.redirect_stdout(io.StringIO()):
        result = instance.enumerate(labels)
    assert instance.active_backend == "python"
    assert len(result) == 257
    assert not instance._use_bytes_key


def test_native_resets_context_and_preserves_domain_input_order(core):
    pools = {3: ["Z", "B", "A"]}
    options = dict(pools=pools, degrees=[3]*4)
    _, expected = run(HashStorageMOFEnumerator, dihedral(4), **options)
    instance, result = run(HybridMOFEnumerator, dihedral(4), backend="cpp", **options)
    with contextlib.redirect_stdout(io.StringIO()):
        repeated = instance.enumerate(["A", "B", "C"])
    assert result == repeated == expected
    assert instance._native_context is None
    assert instance.user_selected_ligands is pools


def test_native_expired_task_does_not_claim_coverage(core, monkeypatch):
    # Expire the root deadline after task setup, before native search starts.
    instance = HybridMOFEnumerator(dihedral(4), 4, backend="cpp", verbose=False, task_time_limit=1e-12)
    with contextlib.redirect_stdout(io.StringIO()):
        result = instance.enumerate(["A", "B"])
    assert not result
    assert instance.task_dedup_stats["completed_tasks"] == 0
    assert not instance._covered_partitions


@pytest.mark.parametrize("kernel", ["dfs", "counter"])
def test_native_overlap_masks_support_more_than_64_bits(core, monkeypatch, kernel):
    _, expected = run(HashStorageMOFEnumerator, dihedral(4))
    original = TaskDedupMOFEnumerator._compile_exclusions
    def shifted(self, *args):
        events, decided, initial = original(self, *args)
        return [[(condition, mask << 70) for condition, mask in level] for level in events], \
            [mask << 70 for mask in decided], initial << 70
    monkeypatch.setattr(TaskDedupMOFEnumerator, "_compile_exclusions", shifted)
    _, result = run(HybridMOFEnumerator, dihedral(4), backend="cpp", native_kernel=kernel)
    assert result == expected


def test_native_mid_search_timeout_retains_partial_results_without_coverage(core):
    n = 20
    instance = HybridMOFEnumerator([tuple(range(n))], n, backend="cpp", verbose=False,
                                   task_time_limit=.01, max_total_results=None,
                                   identity_orbit_limit=float("inf"))
    with contextlib.redirect_stdout(io.StringIO()):
        result = instance.enumerate(["A", "B"])
    assert 0 < len(result) < 2 ** n
    assert any(task["timed_out"] for task in instance.native_task_reports)
    assert not instance._covered_partitions
    assert instance.task_dedup_stats["completed_tasks"] == 0


def test_main_defaults_to_optimized_auto_and_keeps_explicit_legacy():
    from main import get_enumerator_implementation
    from enumerator import MOFEnumerator
    assert get_enumerator_implementation() == (HybridMOFEnumerator, {"backend": "auto"})
    assert get_enumerator_implementation("legacy") == (MOFEnumerator, {})
    for backend in ("python", "auto", "cpp"):
        cls, options = get_enumerator_implementation(backend)
        assert cls is HybridMOFEnumerator
        assert options == {"backend": backend}


def test_cli_allows_native_exact_rho_quota(monkeypatch):
    import main
    class Selected(Exception):
        pass
    def stop_before_input_load(*args, **kwargs):
        raise Selected
    monkeypatch.setattr(main, "load_graph", stop_before_input_load)
    monkeypatch.setattr("sys.argv", ["main.py", "--enumeration-backend", "cpp",
                                     "--rho-quota-per-level", "20", "--prefix-trie"])
    with pytest.raises(Selected):
        main.main()


def test_native_tuple_storage_option_and_zero_limit(core):
    _, expected = run(HashStorageMOFEnumerator, dihedral(4), storage_key="tuple")
    _, actual = run(HybridMOFEnumerator, dihedral(4), backend="cpp", storage_key="tuple")
    assert actual == expected
    _, empty = run(HybridMOFEnumerator, dihedral(4), backend="cpp", max_total_results=0)
    assert empty == []
