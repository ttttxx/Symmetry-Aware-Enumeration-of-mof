"""Compact native result index: retained lookup, byte alphabet and lifecycle."""
import contextlib
import io
import itertools

import pytest

from enumerator_native import HybridMOFEnumerator, native_available
from tests.test_subgroup_task_dedup import dihedral, run


@pytest.fixture(autouse=True)
def core():
    if not native_available():
        pytest.skip("build optional C++ core")


@pytest.mark.parametrize("prefix", [False, True])
@pytest.mark.parametrize("limits", [{}, {"max_total_results": 7},
                                    {"identity_orbit_limit": 0, "identity_new_limit": 2}])
def test_trie_and_hash_keep_ids_records_and_budget_sequence(prefix, limits):
    hashed, expected = run(HybridMOFEnumerator, dihedral(5), backend="cpp",
                           native_store="hash", prefix_trie=prefix, **limits)
    trie, records = run(HybridMOFEnumerator, dihedral(5), backend="cpp",
                        native_store="trie", prefix_trie=prefix, **limits)
    assert records == expected
    assert not trie._canonical_index
    assert trie.result_index_stats["hash_keys"] == 0
    assert trie.result_index_stats["node_bytes"] == 16 * trie.result_index_stats["nodes"]
    assert trie.result_index_stats["nodes"] <= 2 * len(records) + 1
    for record in records:
        assert trie.classes[trie._trie_lookup(tuple(record["coloring"]))] is record
        assert hashed.classes[hashed._trie_lookup(tuple(record["coloring"]))] == record
    assert trie._trie_lookup(("unknown",) * 5) is None
    assert trie._trie_lookup(("A",)) is None
    if limits:
        present = {tuple(r["coloring"]) for r in records}
        missing = next(c for c in itertools.product("ABC", repeat=5) if c not in present)
        assert trie._trie_lookup(missing) is None


def test_unsigned_byte_codes_and_retained_lookup_after_label_decoding():
    labels = [f"label{i:03}" for i in range(256)]
    instance = HybridMOFEnumerator([(0, 1)], 2, backend="cpp", verbose=False,
                                   max_total_results=None, identity_orbit_limit=float("inf"))
    with contextlib.redirect_stdout(io.StringIO()):
        records = instance.enumerate(labels)
    assert len(records) == 256 ** 2
    assert not instance._canonical_index
    for left in (0, 1, 127, 128, 254, 255):
        for right in (0, 1, 127, 128, 254, 255):
            assert instance._trie_lookup((labels[left], labels[right])) == left * 256 + right
    with contextlib.redirect_stdout(io.StringIO()):
        repeated = instance.enumerate(["Y", "Z"])
    assert len(repeated) == 4
    assert instance._trie_lookup(("Y", "Z")) == 1
    assert instance._trie_lookup((labels[0], labels[0])) is None
    assert instance.result_index_stats["nodes"] <= 9


def test_native_storage_selection_is_explicit_in_cli(monkeypatch, capsys):
    import main
    cls, options = main.get_enumerator_implementation("cpp", native_store="trie")
    assert cls is HybridMOFEnumerator
    assert options == dict(backend="cpp", native_store="trie")
    monkeypatch.setattr("sys.argv", ["main.py", "--enumeration-backend", "legacy", "--native-store", "trie"])
    with pytest.raises(SystemExit) as exc:
        main.main()
    assert exc.value.code == 2
    assert "requires --enumeration-backend" in capsys.readouterr().err
    class InputReached(Exception):
        pass
    def stop_before_input(*args, **kwargs):
        raise InputReached
    monkeypatch.setattr(main, "load_graph", stop_before_input)
    monkeypatch.setattr("sys.argv", ["main.py", "--native-store", "trie"])
    with pytest.raises(InputReached):
        main.main()
