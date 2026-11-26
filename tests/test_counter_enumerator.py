"""Counter/DFS equivalence, including preserved subgroup scheduling."""
import contextlib
import io
import itertools

import pytest

from tests.test_subgroup_task_dedup import dihedral, run
from enumerator_casm_inspired import CasmInspiredMOFEnumerator
from enumerator_counter import (CounterMOFEnumerator, SingleTaskCounterMOFEnumerator,
                                SingleTaskDfsMOFEnumerator)


@pytest.mark.parametrize("perms", [dihedral(4), dihedral(5), list(itertools.permutations(range(3)))])
def test_counter_preserves_subgroup_dfs_output_order_and_metadata(perms):
    _, expected = run(CasmInspiredMOFEnumerator, perms)
    instance, actual = run(CounterMOFEnumerator, perms)
    assert actual == expected
    assert instance.task_dedup_stats["duplicate_leaves"] == 0


@pytest.mark.parametrize("perms", [dihedral(4), dihedral(5)])
def test_single_task_counter_and_dfs_have_identical_complete_orbits(perms):
    _, expected = run(CasmInspiredMOFEnumerator, perms)
    _, dfs = run(SingleTaskDfsMOFEnumerator, perms)
    instance, counter = run(SingleTaskCounterMOFEnumerator, perms)
    signature = lambda rows: {tuple(r["coloring"]): (r["stabilizer_size"], r["orbit_size"]) for r in rows}
    assert counter == dfs
    assert signature(counter) == signature(expected)
    assert instance.task_dedup_stats["counter_candidates"] == 3 ** len(perms[0])


@pytest.mark.parametrize("limits", [dict(max_total_results=7),
                                   dict(identity_orbit_limit=0, identity_new_limit=2)])
def test_grouped_counter_keeps_quantity_budget_results(limits):
    _, expected = run(CasmInspiredMOFEnumerator, dihedral(4), **limits)
    _, actual = run(CounterMOFEnumerator, dihedral(4), **limits)
    assert actual == expected


def test_counter_handles_single_choice_domains_without_reordering_labels():
    options = dict(pools={3: ["Z"], 4: ["Y", "X"]}, degrees=[3, 3, 4, 4])
    perms = [(0, 1, 2, 3), (1, 0, 3, 2)]
    _, expected = run(CasmInspiredMOFEnumerator, perms, **options)
    _, actual = run(CounterMOFEnumerator, perms, **options)
    assert actual == expected


def test_counter_does_not_invoke_recursive_dfs_engine(monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("Recursive DFS must not run in the counter experiment")
    monkeypatch.setattr(CasmInspiredMOFEnumerator, "_search_compiled_task", fail)
    instance, result = run(SingleTaskCounterMOFEnumerator, dihedral(4))
    assert result
    assert instance.task_dedup_stats["counter_candidates"] == 81
