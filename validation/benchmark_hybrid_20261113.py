"""Compare Python/C++ DFS and C++ counter with identical subgroup/output work."""
import contextlib
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import resource
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
from enumerator_native import HybridMOFEnumerator


def run_case(topology, c3, c4):
    graph = enrich_graph_data_with_topology(load_graph(str(ROOT / "data" / topology / f"{topology}_graph.json")))
    pools = build_ligand_dict(c3, c4)
    labels = sorted({label for pool in pools.values() for label in pool})
    perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
    kwargs = dict(verbose=False, topology_requirements=calculate_topology_requirements(graph),
                  user_selected_ligands=pools, node_degrees=[node["degree"] for node in graph["nodes"]],
                  max_total_results=None, identity_orbit_limit=float("inf"))
    variants = [("python_dfs", "python", "dfs"), ("cpp_dfs", "cpp", "dfs"),
                ("cpp_counter", "cpp", "counter")]
    times = {name: [] for name, _, _ in variants}
    reports = {}
    expected = None
    digest = None
    for repetition in range(3):
        ordered = variants[repetition:] + variants[:repetition]
        for name, backend, kernel in ordered:
            instance = HybridMOFEnumerator(perms, graph["n_nodes"], backend=backend,
                                           native_kernel=kernel, **kwargs)
            gc.collect()
            started = time.perf_counter()
            results = instance.enumerate(labels)
            elapsed = time.perf_counter()-started
            assert instance.active_backend == backend
            if expected is None:
                expected = results
                digest = hashlib.sha256(json.dumps(results, separators=(",", ":")).encode()).hexdigest()
            assert results == expected
            times[name].append(elapsed)
            reports[name] = dict(native_kernel_s=sum(r["kernel_s"] for r in instance.native_task_reports),
                                  native_tasks=instance.native_task_reports,
                                  task_stats=instance.task_dedup_stats)
    naive = math.prod(len(pools[d]) ** n for d, n in kwargs["topology_requirements"].items())
    assert sum(r["orbit_size"] for r in expected) == naive
    bounded = []
    for limits in [dict(max_total_results=7), dict(max_total_results=100),
                   dict(identity_orbit_limit=0, identity_new_limit=3)]:
        bounded_expected = None
        for name, backend, kernel in variants:
            instance = HybridMOFEnumerator(perms, graph["n_nodes"], backend=backend,
                                           native_kernel=kernel, **dict(kwargs, **limits))
            result = instance.enumerate(labels)
            if bounded_expected is None:
                bounded_expected = result
            assert result == bounded_expected
        bounded.append(dict(limits=limits, count=len(bounded_expected), all_variants_match=True))
    return dict(case=[topology,c3,c4], results=len(expected), output_sha256=digest,
                exact_output_order_and_metadata_match=True, timings_s=times,
                medians_s={name: statistics.median(values) for name, values in times.items()},
                reports=reports, bounded_checks=bounded)


if __name__ == "__main__":
    with contextlib.redirect_stdout(sys.stderr):
        rows = []
        for case in [("ffc",2,2),("ffc",3,3),("fjh",1,2),("fjh",2,1),("fjh",1,3)]:
            row = run_case(*case)
            rows.append(row)
            print("Finished", case, row["medians_s"], file=sys.stderr)
    print(json.dumps(rows, indent=2))
