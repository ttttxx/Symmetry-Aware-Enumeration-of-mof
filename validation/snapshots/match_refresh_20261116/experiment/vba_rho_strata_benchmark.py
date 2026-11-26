#!/usr/bin/env python3
"""Compare capped legacy vba output with exact-rho-layer output quotas.

This enumeration-only benchmark does not assemble CIFs or overwrite the
production result directories. It writes a dated protocol, run logs, compressed
representative lists, and a per-stratum summary under results/.
"""

import csv
import gzip
import json
import os
import platform
import sys
import time
from collections import Counter
from contextlib import redirect_stdout
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from enumerator import MOFEnumerator
from enumerator_rho_quota import RhoQuotaMOFEnumerator
from main import calculate_topology_requirements, enrich_graph_data_with_topology, load_graph
from symmetry import SymmetryEngine


TOPOLOGY = "vba"
OUTPUT_CAP = 100_000
TASK_TIME_LIMIT_S = 60
RHO_LEVELS = (16, 8, 4, 2, 1)
PER_LEVEL_QUOTA = OUTPUT_CAP // len(RHO_LEVELS)
POOLS = {
    "C4x2_Cu_B": {4: ["4c_Cu_1_Ch", "4c_B_Ch"]},
    "C4x3_Cu_B_Cd": {4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"]},
}
OUTPUT_DIR = PROJECT_ROOT / "results" / "vba_rho_stratification_20261010"


def _run(pool_name, ligand_pool, strategy, graph_data, permutations,
         output_cap=OUTPUT_CAP, run_label=None):
    selected_ligands = {degree: list(values) for degree, values in ligand_pool.items()}
    atom_types = sorted({atom for values in selected_ligands.values() for atom in values})
    topology_requirements = calculate_topology_requirements(graph_data)
    node_degrees = [node.get("degree", 0) for node in graph_data["nodes"]]
    options = dict(
        verbose=True,
        topology_requirements=topology_requirements,
        user_selected_ligands=selected_ligands,
        node_degrees=node_degrees,
        task_time_limit=TASK_TIME_LIMIT_S,
        max_total_results=output_cap,
    )
    if strategy == "rho_quota":
        enumerator = RhoQuotaMOFEnumerator(
            permutations, graph_data["n_nodes"],
            rho_quota_per_level=PER_LEVEL_QUOTA,
            **options,
        )
    else:
        enumerator = MOFEnumerator(permutations, graph_data["n_nodes"], **options)

    prefix = f"{pool_name}_{run_label or strategy}"
    log_path = OUTPUT_DIR / f"{prefix}.log"
    t0 = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log, redirect_stdout(log):
        results = enumerator.enumerate(atom_types)
    elapsed_s = time.perf_counter() - t0

    by_order = {}
    for order in RHO_LEVELS:
        members = [item for item in results if item["stabilizer_size"] == order]
        composition_vectors = {
            tuple(Counter(item["coloring"])[atom] for atom in atom_types)
            for item in members
        }
        layer_report = None
        if strategy == "rho_quota":
            layer_report = next(
                (layer for layer in enumerator.rho_layer_report
                 if layer["stabilizer_order"] == order),
                None,
            )
        by_order[order] = {
            "count": len(members),
            "distinct_composition_vectors": len(composition_vectors),
            "layer_report": layer_report,
        }

    total_assignments = 1
    for degree, count in topology_requirements.items():
        total_assignments *= len(selected_ligands.get(degree, [])) ** count
    burnside_rhs = total_assignments / len(permutations)
    burnside_lhs = sum(1.0 / item["stabilizer_size"] for item in results)
    record = {
        "topology": TOPOLOGY,
        "pool": pool_name,
        "ligands_by_degree": selected_ligands,
        "strategy": strategy,
        "output_cap": output_cap,
        "per_level_quota": PER_LEVEL_QUOTA if strategy == "rho_quota" else None,
        "task_time_limit_s": TASK_TIME_LIMIT_S,
        "group_order": len(permutations),
        "returned": len(results),
        "elapsed_s": elapsed_s,
        "rho_setup_s": getattr(enumerator, "rho_setup_time_s", None),
        "burnside_coverage_percent": 100.0 * burnside_lhs / burnside_rhs,
        "rho_layers": by_order,
        "layer_report": getattr(enumerator, "rho_layer_report", None),
        "log": log_path.name,
    }

    raw_path = OUTPUT_DIR / f"{prefix}_representatives.jsonl.gz"
    with gzip.open(raw_path, "wt", encoding="utf-8") as stream:
        for item in results:
            stream.write(json.dumps({
                "coloring": item["coloring"],
                "stabilizer_size": item["stabilizer_size"],
                "orbit_size": item["orbit_size"],
            }, separators=(",", ":")) + "\n")
    record["representatives"] = raw_path.name
    return record


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    graph_path = PROJECT_ROOT / "data" / TOPOLOGY / f"{TOPOLOGY}_graph.json"
    graph_data = enrich_graph_data_with_topology(load_graph(str(graph_path)))
    permutations = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()

    protocol = {
        "purpose": "Compare output allocation across exact stabilizer-order strata under the existing 100,000-result cap; not a sample-quality or symmetry-preference test.",
        "topology": TOPOLOGY,
        "pools": POOLS,
        "group_order": len(permutations),
        "rho_orders": list(RHO_LEVELS),
        "output_cap": OUTPUT_CAP,
        "per_level_quota": PER_LEVEL_QUOTA,
        "task_time_limit_s": TASK_TIME_LIMIT_S,
        "python": sys.version,
        "platform": platform.platform(),
        "methods": [
            "legacy_subgroup_dfs_at_100k_cap",
            "exact_rho_layer_quota_at_100k_cap",
            "legacy_subgroup_dfs_matched_to_rho_quota_return_count",
        ],
        "measurement": "enumeration wall time only; graph loading and automorphism discovery excluded",
    }
    (OUTPUT_DIR / "protocol.json").write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    records = []
    for pool_name, ligand_pool in POOLS.items():
        for strategy in ("legacy", "rho_quota"):
            print(f"Starting {pool_name} / {strategy}; outputs -> {OUTPUT_DIR}", flush=True)
            record = _run(
                pool_name, ligand_pool, strategy, graph_data, permutations,
                output_cap=OUTPUT_CAP,
            )
            records.append(record)
            (OUTPUT_DIR / f"{pool_name}_{strategy}_summary.json").write_text(
                json.dumps(record, indent=2), encoding="utf-8"
            )
            print(
                f"Finished {pool_name} / {strategy}: {record['returned']} outputs "
                f"in {record['elapsed_s']:.1f}s", flush=True
            )
            if strategy == "rho_quota":
                matched_cap = record["returned"]
                print(
                    f"Starting {pool_name} / legacy_matched_{matched_cap}; "
                    f"outputs -> {OUTPUT_DIR}", flush=True,
                )
                matched = _run(
                    pool_name, ligand_pool, "legacy", graph_data, permutations,
                    output_cap=matched_cap,
                    run_label=f"legacy_matched_{matched_cap}",
                )
                records.append(matched)
                (OUTPUT_DIR / f"{pool_name}_legacy_matched_summary.json").write_text(
                    json.dumps(matched, indent=2), encoding="utf-8"
                )
                print(
                    f"Finished matched legacy {pool_name}: {matched['returned']} "
                    f"outputs in {matched['elapsed_s']:.1f}s", flush=True,
                )

    csv_path = OUTPUT_DIR / "rho_layer_summary.csv"
    fields = [
        "pool", "strategy", "output_cap", "per_level_quota", "returned",
        "elapsed_s", "rho_order", "rho", "layer_count",
        "distinct_composition_vectors", "layer_elapsed_s", "layer_status",
        "tasks_started", "fixed_layer_tasks", "task_timeouts",
        "burnside_coverage_percent",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for record in records:
            for order, layer in record["rho_layers"].items():
                report = layer["layer_report"] or {}
                writer.writerow({
                    "pool": record["pool"],
                    "strategy": record["strategy"],
                    "output_cap": record["output_cap"],
                    "per_level_quota": record["per_level_quota"],
                    "returned": record["returned"],
                    "elapsed_s": round(record["elapsed_s"], 3),
                    "rho_order": order,
                    "rho": f"{int(order) / record['group_order']:.6g}",
                    "layer_count": layer["count"],
                    "distinct_composition_vectors": layer["distinct_composition_vectors"],
                    "layer_elapsed_s": round(report.get("elapsed_s", 0), 3)
                    if report else "",
                    "layer_status": report.get("status", "not_layered"),
                    "tasks_started": report.get("tasks_started", ""),
                    "fixed_layer_tasks": report.get("fixed_layer_tasks", ""),
                    "task_timeouts": report.get("task_timeouts", ""),
                    "burnside_coverage_percent": round(record["burnside_coverage_percent"], 6),
                })
    print(f"Summary written to {csv_path}", flush=True)


if __name__ == "__main__":
    main()
