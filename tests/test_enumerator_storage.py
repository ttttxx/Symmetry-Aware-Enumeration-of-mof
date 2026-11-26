"""Storage replacements must preserve the full subgroup DFS output contract."""
import contextlib
import io

import pytest

from tests.test_subgroup_task_dedup import dihedral, run
from enumerator_casm_inspired import CasmInspiredMOFEnumerator
from enumerator_storage import HashStorageMOFEnumerator


@pytest.mark.parametrize("key", ["tuple", "bytes"])
@pytest.mark.parametrize("limits", [dict(), dict(max_total_results=7),
                                   dict(identity_orbit_limit=0, identity_new_limit=2)])
def test_storage_preserves_complete_and_bounded_output(key, limits):
    _, expected = run(CasmInspiredMOFEnumerator, dihedral(4), **limits)
    instance, result = run(HashStorageMOFEnumerator, dihedral(4), storage_key=key, **limits)
    assert result == expected
    assert len(instance._canonical_index) == len(result)
    # IDs in the membership index refer to insertion positions, as in the Trie.
    assert instance._trie_lookup(tuple(instance.classes[0]["coloring"])) == 0
    assert instance._trie_lookup(("unknown",) * 4) is None


@pytest.mark.parametrize("key", ["tuple", "bytes"])
def test_storage_resets_between_calls_and_keeps_original_label_order(key):
    pools = {3: ["Z", "A", "B"]}
    options = dict(pools=pools, degrees=[3] * 4)
    _, expected = run(CasmInspiredMOFEnumerator, dihedral(4), **options)
    instance, first = run(HashStorageMOFEnumerator, dihedral(4), storage_key=key, **options)
    with contextlib.redirect_stdout(io.StringIO()):
        second = instance.enumerate(["A", "B", "C"])
    assert first == second == expected
    assert instance.user_selected_ligands is pools


def test_byte_key_falls_back_to_exact_tuple_above_256_labels():
    # One site keeps the independent exhaustive check small despite many types.
    labels = [f"type_{i:03}" for i in range(257)]
    kwargs = dict(verbose=False, max_total_results=None, identity_orbit_limit=float("inf"))
    instance = HashStorageMOFEnumerator([(0,)], 1, **kwargs)
    baseline = CasmInspiredMOFEnumerator([(0,)], 1, **kwargs)
    with contextlib.redirect_stdout(io.StringIO()):
        expected = baseline.enumerate(labels)
        actual = instance.enumerate(labels)
    assert actual == expected
    assert not instance._use_bytes_key
    assert all(isinstance(key, tuple) for key in instance._canonical_index)


def test_storage_returns_existing_zero_id_for_duplicate_leaf():
    instance, result = run(HashStorageMOFEnumerator, dihedral(4))
    record = instance.classes[0]
    coloring = tuple(record["coloring"])
    size = len(instance.classes)
    instance._store_canonical_leaf(coloring, record["stabilizer_size"])
    assert len(instance.classes) == size
    assert instance.task_dedup_stats["duplicate_leaves"] == 1
