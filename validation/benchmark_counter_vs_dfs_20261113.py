"""2x2 experiment: full-product counter/DFS x subgroup/single-task schedule."""
import contextlib
import json
import math
import os
from pathlib import Path
import statistics
import sys
import time

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import sympy.combinatorics
from experiment._ours_worker import build_ligand_dict
from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
from symmetry import SymmetryEngine
from enumerator_casm_inspired import CasmInspiredMOFEnumerator
from enumerator_counter import CounterMOFEnumerator, SingleTaskCounterMOFEnumerator, SingleTaskDfsMOFEnumerator


def run_case(topology, c3, c4):
    graph = enrich_graph_data_with_topology(load_graph(str(ROOT / "data" / topology / f"{topology}_graph.json")))
    perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
    pools = build_ligand_dict(c3, c4)
    labels = sorted({label for pool in pools.values() for label in pool})
    kwargs = dict(verbose=False, topology_requirements=calculate_topology_requirements(graph),
                  user_selected_ligands=pools, node_degrees=[node["degree"] for node in graph["nodes"]],
                  max_total_results=None, identity_orbit_limit=float("inf"))
    variants = [("subgroup_dfs", CasmInspiredMOFEnumerator), ("subgroup_counter", CounterMOFEnumerator),
                ("single_dfs", SingleTaskDfsMOFEnumerator), ("single_counter", SingleTaskCounterMOFEnumerator)]
    timings = {name: [] for name, _ in variants}
    stats = {}
    output_by_variant = {}
    signature_expected = None
    # Four repetitions rotate every version through every run position once.
    for repetition in range(4):
        ordered = variants[repetition:] + variants[:repetition]
        for name, cls in ordered:
            instance = cls(perms, graph["n_nodes"], **kwargs)
            started = time.perf_counter()
            result = instance.enumerate(labels)
            elapsed = time.perf_counter()-started
            signature = {tuple(r["coloring"]): (r["stabilizer_size"], r["orbit_size"]) for r in result}
            if signature_expected is None:
                signature_expected = signature
            assert signature == signature_expected
            if name in output_by_variant:
                assert result == output_by_variant[name]
            else:
                output_by_variant[name] = result
            timings[name].append(elapsed)
            stats[name] = instance.task_dedup_stats
        assert output_by_variant["subgroup_counter"] == output_by_variant["subgroup_dfs"]
        assert output_by_variant["single_counter"] == output_by_variant["single_dfs"]
    naive = math.prod(len(pools[d]) ** count for d, count in kwargs["topology_requirements"].items())
    assert sum(r["orbit_size"] for r in output_by_variant["single_counter"]) == naive
    assert stats["single_counter"]["counter_candidates"] == naive
    return dict(case=[topology, c3, c4], naive=naive, representatives=len(signature_expected),
                timings_s=timings, medians_s={name: statistics.median(values) for name, values in timings.items()},
                stats=stats, all_representatives_and_metadata_match=True,
                traversal_order_matches_within_each_schedule=True)


if __name__ == "__main__":
    with contextlib.redirect_stdout(sys.stderr):
        results = []
        for case in [("ffc", 2, 2), ("ffc", 3, 3), ("fjh", 1, 2), ("fjh", 2, 1), ("fjh", 1, 3)]:
            row = run_case(*case)
            results.append(row)
            print("Finished", case, row["medians_s"], file=sys.stderr)
    print(json.dumps(results, indent=2))
