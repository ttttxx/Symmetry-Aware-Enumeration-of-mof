"""Coverage, output order and truncation regressions for the overlap prototype."""
import contextlib
import io
import itertools
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from enumerator import MOFEnumerator
from enumerator_task_dedup import TaskDedupMOFEnumerator
from enumerator_casm_inspired import CasmInspiredMOFEnumerator


def dihedral(n):
    return sorted({tuple((sign * i + shift) % n for i in range(n))
                   for sign in (-1, 1) for shift in range(n)})


def run(cls, perms, pools=None, degrees=None, **limits):
    instance = cls(perms, len(perms[0]), verbose=False,
                   user_selected_ligands=pools, node_degrees=degrees,
                   max_total_results=limits.pop("max_total_results", None),
                   identity_orbit_limit=limits.pop("identity_orbit_limit", float("inf")),
                   **limits)
    with contextlib.redirect_stdout(io.StringIO()):
        result = instance.enumerate(["A", "B", "C"])
    return instance, result


@pytest.mark.parametrize("perms", [dihedral(4), dihedral(5), list(itertools.permutations(range(3)))])
@pytest.mark.parametrize("implementation", [TaskDedupMOFEnumerator, CasmInspiredMOFEnumerator])
def test_overlapping_non_nested_subgroups_preserve_all_orbits(perms, implementation):
    old, expected = run(MOFEnumerator, perms)
    new, actual = run(implementation, perms)
    assert actual == expected  # Includes stable output order and metadata.
    # Independent exhaustive orbit oracle, not just agreement with baseline.
    representatives = {}
    for coloring in itertools.product("ABC", repeat=len(perms[0])):
        canonical = min(tuple(coloring[p.index(i)] for i in range(len(p))) for p in perms)
        representatives[canonical] = sum(tuple(coloring[p.index(i)] for i in range(len(p))) == coloring for p in perms)
    assert {tuple(r["coloring"]): r["stabilizer_size"] for r in actual} == representatives
    assert sum(r["orbit_size"] for r in actual) == 3 ** len(perms[0])
    assert new.task_dedup_stats["leaf_visits"] == len(actual)
    assert new.task_dedup_stats["duplicate_leaves"] == 0
    assert new.task_dedup_stats["pruned_prefixes"] > 0


@pytest.mark.parametrize("implementation", [TaskDedupMOFEnumerator, CasmInspiredMOFEnumerator])
def test_singleton_domains_keep_distinct_colors_and_variable_choices(implementation):
    perms = [(0, 1, 2, 3), (1, 0, 3, 2)]
    options = dict(pools={3: ["A"], 4: ["B", "C"]}, degrees=[3, 3, 4, 4])
    _, expected = run(MOFEnumerator, perms, **options)
    new, actual = run(implementation, perms, **options)
    assert actual == expected
    assert new.task_dedup_stats["leaf_visits"] == len(actual)


@pytest.mark.parametrize("limits", [dict(max_total_results=7),
                                   dict(identity_orbit_limit=0, identity_new_limit=2)])
@pytest.mark.parametrize("implementation", [TaskDedupMOFEnumerator, CasmInspiredMOFEnumerator])
def test_output_cap_and_identity_budget_preserve_sequence(limits, implementation):
    _, expected = run(MOFEnumerator, dihedral(4), **limits)
    _, actual = run(implementation, dihedral(4), **limits)
    assert actual == expected


@pytest.mark.parametrize("implementation", [TaskDedupMOFEnumerator, CasmInspiredMOFEnumerator])
def test_timed_out_task_does_not_exclude_unvisited_candidates(implementation):
    instance = implementation(dihedral(4), 4, verbose=False, max_total_results=None)
    instance._covered_partitions = []
    instance.task_dedup_stats = dict(tasks=0, completed_tasks=0, pruned_prefixes=0,
                                     leaf_visits=0, duplicate_leaves=0, compile_s=0.)
    full = instance._compute_orbits_by_gen(instance.perms)
    args = (["A", "B"], [None] * 4, list(range(1, len(instance.perms))))
    instance._dfs_lattice_pruning(0, full, *args, deadline=time.monotonic()-1)
    assert not instance._covered_partitions
    instance._dfs_lattice_pruning(0, [[i] for i in range(4)], *args)
    assert len(instance.classes) == 6
    assert instance.task_dedup_stats["completed_tasks"] == 1


@pytest.mark.parametrize("implementation", [TaskDedupMOFEnumerator, CasmInspiredMOFEnumerator])
def test_same_instance_can_be_reused_without_stale_exclusions(implementation):
    instance, first = run(implementation, dihedral(4))
    with contextlib.redirect_stdout(io.StringIO()):
        second = instance.enumerate(["A", "B", "C"])
    assert first == second


@pytest.mark.parametrize("implementation", [TaskDedupMOFEnumerator, CasmInspiredMOFEnumerator])
def test_partial_timeout_keeps_remaining_region_searchable(monkeypatch, implementation):
    instance = implementation(dihedral(4), 4, verbose=False, max_total_results=None)
    instance._covered_partitions = []
    instance.task_dedup_stats = dict(tasks=0, completed_tasks=0, pruned_prefixes=0,
                                     leaf_visits=0, duplicate_leaves=0, compile_s=0.)
    ticks = iter([0., 0., 0.])
    with monkeypatch.context() as patch:
        patch.setattr(time, "monotonic", lambda: next(ticks, 2.))
        instance._dfs_lattice_pruning(0, [[0, 1, 2, 3]], ["A", "B"],
                                     [None] * 4, list(range(1, len(instance.perms))), deadline=1.)
    assert len(instance.classes) == 1
    assert not instance._covered_partitions
    instance._dfs_lattice_pruning(0, [[i] for i in range(4)], ["A", "B"],
                                 [None] * 4, list(range(1, len(instance.perms))))
    assert len(instance.classes) == 6


def test_integer_encoding_preserves_unsorted_domain_order_and_public_labels():
    pools = {3: ["zinc", "boron", "copper"]}
    options = dict(pools=pools, degrees=[3] * 4)
    _, expected = run(MOFEnumerator, dihedral(4), **options)
    instance, actual = run(CasmInspiredMOFEnumerator, dihedral(4), **options)
    assert actual == expected
    assert instance.user_selected_ligands is pools
    assert all(instance._trie_lookup(tuple(r["coloring"])) is not None for r in actual)


def test_lazy_canonical_fallback_matches_full_group_images():
    instance = CasmInspiredMOFEnumerator(dihedral(4), 4, verbose=False)
    for coloring in itertools.product("AB", repeat=4):
        assert instance._get_canonical(coloring) == MOFEnumerator._get_canonical(instance, coloring)
        assert instance._stabilizer_size(coloring) == MOFEnumerator._compute_stabilizer_and_orbit(instance, coloring)[0]
