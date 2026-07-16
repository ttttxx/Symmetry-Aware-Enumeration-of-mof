# !/usr/bin/env python3
"""Our enumerator worker — runs under project python3, calls MOFEnumerator directly.

Usage:
    python3 _ours_worker.py <topology> <c3_count> <c4_count> [--disable-budget]

Prints a JSON record to stdout:
    {"orbits": N, "time_s": t, "peak_mem_mb": m, "budget_triggered": bool, "status": "exhaustive"|"budgeted"}
"""

import json
import os
import sys
import time
import resource

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator

DEGREE_LIGANDS = {
    3: ["3c_B_Ch", "3c_triazine_Ch", "3c_Cu_imadazoline_Ch"],
    4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"],
}

def build_ligand_dict(c3_count: int, c4_count: int) -> dict:
    """Select first |C3| and first |C4| ligands from the available pools."""
    ligands = {}
    if c3_count > 0:
        ligands[3] = DEGREE_LIGANDS[3][:c3_count]
    if c4_count > 0:
        ligands[4] = DEGREE_LIGANDS[4][:c4_count]
    return ligands

def ligand_pool_label(ligands: dict) -> str:
    parts = []
    for d in sorted(ligands):
        parts.append(f"C{d}={len(ligands[d])}")
    return ",".join(parts)

def run_ours_enumeration(topology: str, c3_count: int, c4_count: int,
                         disable_budget: bool = False) -> dict:
    """Run our enumerator in-process, measure timing + RSS around enumerate() only."""
    graph_path = os.path.join(PROJECT_ROOT, "data", topology, f"{topology}_graph.json")
    graph_data = load_graph(graph_path)
    graph_data = enrich_graph_data_with_topology(graph_data)
    node_degrees = [n.get("degree", 0) for n in graph_data["nodes"]]
    topo_req = calculate_topology_requirements(graph_data)

    ligands = build_ligand_dict(c3_count, c4_count)
    all_ligs = sorted(set(l for ls in ligands.values() for l in ls))

    sym = SymmetryEngine(graph_data, verbose=False)
    perms = sym.compute_automorphism_group()

    if disable_budget:
        enum_task_time_limit = None
        enum_max_total = None
        enum_orbit_limit = float("inf")
        enum_new_limit = None
    else:
        enum_task_time_limit = 60       # default from main.py
        enum_max_total = 100_000        # default from main.py
        enum_orbit_limit = None         # use module-constant default
        enum_new_limit = None           # use module-constant default

    enumerator = MOFEnumerator(
        perms, graph_data["n_nodes"], verbose=False,
        topology_requirements=topo_req,
        user_selected_ligands=ligands,
        node_degrees=node_degrees,
        task_time_limit=enum_task_time_limit,
        max_total_results=enum_max_total,
        identity_orbit_limit=enum_orbit_limit,
        identity_new_limit=enum_new_limit,
    )

    t0 = time.perf_counter()
    results = enumerator.enumerate(all_ligs)
    t1 = time.perf_counter()
    rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    # Check if budget was triggered (Burnside coverage < 100%)
    lhs = sum(1.0 / r["stabilizer_size"] for r in results)
    total_schemes = 1
    for deg, cnt in topo_req.items():
        total_schemes *= len(ligands.get(deg, all_ligs)) ** cnt
    rhs = total_schemes / len(perms)
    budget_triggered = abs(lhs - rhs) > 1e-6

    return {
        "orbits": len(results),
        "time_s": round(t1 - t0, 6),
        "peak_mem_mb": round(rss_bytes / 1_048_576, 2),
        "budget_triggered": budget_triggered,
        "status": "budgeted" if budget_triggered else "exhaustive",
    }

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("topology")
    parser.add_argument("c3_count", type=int)
    parser.add_argument("c4_count", type=int)
    parser.add_argument("--disable-budget", action="store_true")
    args = parser.parse_args()

    # Redirect all module logging from stdout → stderr so only JSON goes to stdout
    _real_stdout = sys.stdout
    sys.stdout = sys.stderr

    try:
        result = run_ours_enumeration(
            args.topology, args.c3_count, args.c4_count,
            disable_budget=args.disable_budget,
        )
    finally:
        sys.stdout = _real_stdout

    # Print clean JSON to (the real) stdout
    print(json.dumps(result))

if __name__ == "__main__":
    main()
