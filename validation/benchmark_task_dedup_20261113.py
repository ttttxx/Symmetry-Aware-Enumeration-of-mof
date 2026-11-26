"""Warm, interleaved baseline/prototype timing with full-result equality checks."""
import contextlib
import io
import json
import os
from pathlib import Path
import statistics
import sys
import time

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import sympy.combinatorics  # Explicit warm-import boundary, identical for both.
from experiment._ours_worker import build_ligand_dict
from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator
from enumerator_task_dedup import TaskDedupMOFEnumerator
from enumerator_casm_inspired import CasmInspiredMOFEnumerator


def run_case(topology, c3, c4, include_casm_inspired=False, bounded_check=False):
    graph = enrich_graph_data_with_topology(load_graph(str(ROOT / "data" / topology / f"{topology}_graph.json")))
    perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
    pools = build_ligand_dict(c3, c4)
    labels = sorted({label for values in pools.values() for label in values})
    kwargs = dict(verbose=False, topology_requirements=calculate_topology_requirements(graph),
                  user_selected_ligands=pools, node_degrees=[node["degree"] for node in graph["nodes"]],
                  max_total_results=None, identity_orbit_limit=float("inf"))
    if bounded_check:
        checks = []
        for limits in [dict(max_total_results=7), dict(max_total_results=100),
                       dict(identity_orbit_limit=0, identity_new_limit=3)]:
            bounded_kwargs = dict(kwargs, **limits)
            baseline = MOFEnumerator(perms, graph["n_nodes"], **bounded_kwargs).enumerate(labels)
            enhanced = CasmInspiredMOFEnumerator(perms, graph["n_nodes"], **bounded_kwargs)
            result = enhanced.enumerate(labels)
            assert result == baseline, "Bounded result order or metadata changed"
            checks.append(dict(limits=limits, returned=len(result), exact_match=True,
                               task_stats=enhanced.task_dedup_stats))
        return dict(case=[topology, c3, c4], bounded_checks=checks)
    variant_specs = [("baseline", MOFEnumerator), ("task_dedup", TaskDedupMOFEnumerator)]
    if include_casm_inspired:
        variant_specs.append(("casm_inspired", CasmInspiredMOFEnumerator))
    records = {name: [] for name, _ in variant_specs}
    diagnostic_stats = {}
    expected = None
    for repetition in range(3):
        # Rotate execution order so each variant occupies each position once.
        offset = repetition % len(variant_specs)
        variants = variant_specs[offset:] + variant_specs[:offset]
        for name, cls in variants:
            instance = cls(perms, graph["n_nodes"], **kwargs)
            started = time.perf_counter()
            result = instance.enumerate(labels)
            elapsed = time.perf_counter()-started
            if expected is None:
                expected = result
            assert result == expected, "Representatives, metadata or output order changed"
            records[name].append(elapsed)
            if name != "baseline":
                stats = instance.task_dedup_stats
                assert stats["leaf_visits"] == len(result)
                assert stats["duplicate_leaves"] == 0
                diagnostic_stats[name] = stats
    naive = __import__("math").prod(len(pools[d]) ** n for d, n in kwargs["topology_requirements"].items())
    assert sum(r["orbit_size"] for r in expected) == naive
    return dict(case=[topology, c3, c4], representatives=len(expected), timings_s=records,
                medians_s={name: statistics.median(values) for name, values in records.items()},
                task_dedup_stats=diagnostic_stats["task_dedup"],
                variant_stats=diagnostic_stats,
                exact_order_and_metadata_match=True, covered_assignments=naive)


if __name__ == "__main__":
    with contextlib.redirect_stdout(sys.stderr):
        results = []
        for case in [("ffc", 2, 2), ("ffc", 3, 3), ("fjh", 1, 2), ("fjh", 2, 1), ("fjh", 1, 3)]:
            results.append(run_case(*case, include_casm_inspired="--casm-inspired" in sys.argv,
                                    bounded_check="--bounded-check" in sys.argv))
            print("Finished", case, results[-1].get("medians_s", "bounded checks matched"), file=sys.stderr)
    print(json.dumps(results, indent=2))
