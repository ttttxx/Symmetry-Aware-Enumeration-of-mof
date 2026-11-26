"""Sweep X-site fusion tolerances over the existing 27-structure audit set."""

import argparse
import csv
import os
import random
import sys
import numpy as np
from pymatgen.core import Structure

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from assembler import StructureAssembler
from experiment.geometric_reconstruction_summary import (
    CONFIGS, _enumerate_case, _sample_case, RANDOM_SEED,
)
from structure_validator import StructureValidator

THRESHOLDS = (0.4, 0.6, 0.8, 1.0, 1.2)
SOURCE_RECORDS = os.path.join(
    PROJECT_ROOT, "results", "reconstruction_summary_split",
    "partial_2026-08-28", "reconstruction_records.csv",
)


def load_audit_cases():
    with open(SOURCE_RECORDS, newline="", encoding="utf-8") as stream:
        saved_rows = list(csv.DictReader(stream))
    cases = []
    for config_idx, config in enumerate(CONFIGS[:2]):
        topology = config["topology"]
        graph_data, perms, classified = _enumerate_case(config)
        sampled = _sample_case(
            classified, gamma_size=len(perms),
            sampling=config.get("sampling", {}), seed=RANDOM_SEED + config_idx,
        )
        by_id = {int(item["result"]["id"]): item for item in sampled}
        for row in saved_rows:
            if row["topology"] != topology:
                continue
            structure_id = int(row["structure_id"])
            if structure_id not in by_id:
                raise ValueError(f"Cannot recover {topology} structure {structure_id}")
            rec = by_id[structure_id]
            if rec["sym_class"] != row["sym_class"] or abs(rec["rho"] - float(row["rho"])) > 1e-12:
                raise ValueError(f"Audit metadata mismatch for {topology} structure {structure_id}")
            cases.append((config, graph_data, rec, row))
    if len(cases) != 27:
        raise ValueError(f"Expected 27 audit cases, found {len(cases)}")
    return cases


def run(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    cases = load_audit_cases()
    rows = []
    for tolerance in THRESHOLDS:
        print(f"Threshold {tolerance:.1f} Å: reconstructing {len(cases)} audit cases", flush=True)
        for index, (config, graph_data, rec, source_row) in enumerate(cases):
            topology = config["topology"]
            config_index = next(i for i, item in enumerate(CONFIGS[:2]) if item["topology"] == topology)
            placement_seed = RANDOM_SEED + config_index * 100000 + int(source_row["structure_id"])
            random.seed(placement_seed)
            np.random.seed(placement_seed)
            out = os.path.join(output_dir, f"{topology}_{source_row['structure_id']}_{tolerance:.1f}.cif")
            assembler = StructureAssembler(
                output_dir=output_dir,
                validator=StructureValidator(strict_mode=False, check_symmetry=False),
                run_symmetry=False,
                x_site_fusion_tolerance=tolerance,
            )
            row = {
                "topology": topology,
                "structure_id": source_row["structure_id"],
                "tolerance_A": tolerance,
                "source_cif": source_row["output_cif"],
                "assembly_ok": False,
                "topology_linked_pairs_merged": None,
                "unlinked_pairs_merged": None,
                "topology_linked_pairs": None,
                "merged_pairs": None,
                "topology_pass": None,
                "clash_pass": None,
                "cif_reload_ok": False,
                "cif_n_sites": None,
                "error": "",
            }
            try:
                coloring = rec["result"]["coloring"]
                metadata = assembler.assemble_and_save(
                    filename=out,
                    structure_map={str(i): color for i, color in enumerate(coloring)},
                    graph_data=graph_data,
                    forced_edge=config["forced_edge"],
                    return_metadata=True,
                )
                audit = metadata["x_site_fusion"]
                merged_pairs = set(tuple(pair) for pair in audit["merged_pair_ids"])
                linked_pairs = set(tuple(pair) for pair in audit["topology_linked_pair_ids"])
                row.update({
                    "assembly_ok": True,
                    "topology_linked_pairs_merged": len(merged_pairs & linked_pairs),
                    "unlinked_pairs_merged": len(merged_pairs - linked_pairs),
                    "topology_linked_pairs": len(linked_pairs),
                    "merged_pairs": audit["merged_pairs"],
                    "topology_pass": metadata["validation"].get("topology", {}).get("passed"),
                    "clash_pass": metadata["validation"].get("clash", {}).get("passed"),
                })
                loaded = Structure.from_file(out)
                row["cif_reload_ok"] = True
                row["cif_n_sites"] = len(loaded)
            except Exception as error:
                row["error"] = f"{type(error).__name__}: {error}"
            rows.append(row)
            if (index + 1) % 5 == 0:
                print(f"  {index + 1}/{len(cases)} completed", flush=True)
    csv_path = os.path.join(output_dir, "x_site_fusion_sensitivity.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path = os.path.join(output_dir, "x_site_fusion_sensitivity_summary.csv")
    summaries = []
    for tolerance in THRESHOLDS:
        subset = [r for r in rows if r["tolerance_A"] == tolerance]
        summaries.append({
            "tolerance_A": tolerance,
            "samples": len(subset),
            "assembled": sum(bool(r["assembly_ok"]) for r in subset),
            "topology_linked_pairs_merged": sum(r["topology_linked_pairs_merged"] or 0 for r in subset),
            "unlinked_pairs_merged": sum(r["unlinked_pairs_merged"] or 0 for r in subset),
            "topology_pass": sum(r["topology_pass"] is True for r in subset),
            "cif_reload_ok": sum(bool(r["cif_reload_ok"]) for r in subset),
        })
    with open(summary_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    return csv_path, summary_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=os.path.join(
        PROJECT_ROOT, "results", "x_site_fusion_sensitivity_2026-10-10"))
    args = parser.parse_args()
    for path in run(args.output_dir):
        print(path)
