#!/usr/bin/env python3
"""Trace a small ffc CO2/N2 assembler repair regression set."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, ".")

from assembler import StructureAssembler
from main import enrich_graph_data_with_topology, load_graph


DEFAULT_MANIFEST = Path("results/co2n2_ffc_screening/batch_001_with_validation/manifest.csv")
DEFAULT_OUT_DIR = Path("results/co2n2_ffc_screening/ffc_repair_regression_trace")
FORCED_EDGE = "1B_4H_Ch"

REGRESSION_IDS = [
    # clean
    "ffc_co2n2_batch001_0142",
    "ffc_co2n2_batch001_0143",
    "ffc_co2n2_batch001_0152",
    # clash-only
    "ffc_co2n2_batch001_0173",
    "ffc_co2n2_batch001_0056",
    "ffc_co2n2_batch001_0057",
    # topology-only
    "ffc_co2n2_batch001_0147",
    # clash + topology, covering severe 3/8, 4/8, 5/8, and 6/8 recovery cases
    "ffc_co2n2_batch001_0170",
    "ffc_co2n2_batch001_0132",
    "ffc_co2n2_batch001_0116",
    "ffc_co2n2_batch001_0005",
]


CSV_FIELDS = [
    "candidate_id",
    "source_enum_id",
    "baseline_warning",
    "new_warning",
    "coloring_json",
    "scale_factor",
    "lattice_scaled",
    "merge_count",
    "bridge_bond_count",
    "bridge_repair_count",
    "bridge_repair_stats",
    "edge_rotation_count",
    "edge_rotation_trace",
    "connection_closure",
    "n_clashes",
    "worst_clash_ratio",
    "clash_message",
    "clash_source_summary",
    "clash_source_trace",
    "topology_n_nodes",
    "topology_n_edges",
    "topology_message",
    "atom_trace",
    "prealign_trace",
    "cif_path",
]


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _read_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return {row["candidate_id"]: row for row in csv.DictReader(handle)}


def _failed_validation_checks(validation: dict[str, Any]) -> list[str]:
    failed = []
    for name in ("clash", "topology", "symmetry"):
        result = validation.get(name, {})
        if isinstance(result, dict) and not result.get("skipped") and result.get("passed") is False:
            failed.append(name)
    return failed


def _summarize_clash_sources(trace: dict[str, Any]) -> str:
    summary = trace.get("summary", {}) if isinstance(trace, dict) else {}
    counts = summary.get("category_counts", {}) if isinstance(summary, dict) else {}
    if not isinstance(counts, dict) or not counts:
        return ""
    return "; ".join(f"{key}={value}" for key, value in sorted(counts.items()))


def _trace_row(row: dict[str, str], metadata: dict[str, Any], cif_path: Path) -> dict[str, Any]:
    validation = metadata.get("validation", {})
    clash = validation.get("clash", {}) if isinstance(validation, dict) else {}
    topology = validation.get("topology", {}) if isinstance(validation, dict) else {}
    worst = clash.get("worst_clash", {}) if isinstance(clash, dict) else {}
    if not isinstance(worst, dict):
        worst = {}
    clash_source_trace = metadata.get("clash_source_trace", {})
    new_warning = ";".join(_failed_validation_checks(validation)) if isinstance(validation, dict) else ""
    return {
        "candidate_id": row["candidate_id"],
        "source_enum_id": row["source_enum_id"],
        "baseline_warning": row.get("validation_warning", ""),
        "new_warning": new_warning,
        "coloring_json": row["coloring_json"],
        "scale_factor": metadata.get("scale_factor"),
        "lattice_scaled": metadata.get("lattice_scaled"),
        "merge_count": metadata.get("merge_count"),
        "bridge_bond_count": metadata.get("bridge_bond_count"),
        "bridge_repair_count": metadata.get("bridge_repair_count"),
        "bridge_repair_stats": metadata.get("bridge_repair_stats"),
        "edge_rotation_count": len(metadata.get("edge_rotation_trace") or []),
        "edge_rotation_trace": metadata.get("edge_rotation_trace"),
        "connection_closure": metadata.get("connection_closure"),
        "n_clashes": clash.get("n_clashes") if isinstance(clash, dict) else None,
        "worst_clash_ratio": worst.get("ratio"),
        "clash_message": clash.get("message") if isinstance(clash, dict) else "",
        "clash_source_summary": _summarize_clash_sources(clash_source_trace),
        "clash_source_trace": clash_source_trace,
        "topology_n_nodes": f"{topology.get('n_assembled')}/{topology.get('n_tg')}"
        if isinstance(topology, dict) else "",
        "topology_n_edges": f"{topology.get('e_assembled')}/{topology.get('e_tg')}"
        if isinstance(topology, dict) else "",
        "topology_message": topology.get("message") if isinstance(topology, dict) else "",
        "atom_trace": metadata.get("atom_trace"),
        "prealign_trace": metadata.get("prealign_trace"),
        "cif_path": str(cif_path),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                field: json.dumps(row[field], ensure_ascii=False, sort_keys=True, default=_json_default)
                if isinstance(row.get(field), (dict, list)) else row.get(field, "")
                for field in CSV_FIELDS
            })


def _write_markdown(path: Path, rows: list[dict[str, Any]]) -> None:
    lines = [
        "# ffc Repair Regression Trace",
        "",
        "| candidate | enum | baseline | new | scale | bridges | closure rmsd | clashes | worst ratio | clash sources | topology edges |",
        "|---|---:|---|---|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        closure = row.get("connection_closure") or {}
        closure_rmsd = closure.get("closure_rmsd_A") if isinstance(closure, dict) else None
        lines.append(
            "| {candidate_id} | {source_enum_id} | {baseline_warning} | {new_warning} | "
            "{scale_factor:.4f} | {bridge_bond_count} | {closure} | {n_clashes} | "
            "{worst} | {clash_sources} | {topology_n_edges} |".format(
                candidate_id=row["candidate_id"],
                source_enum_id=row["source_enum_id"],
                baseline_warning=row["baseline_warning"] or "clean",
                new_warning=row["new_warning"] or "clean",
                scale_factor=float(row["scale_factor"] or 0.0),
                bridge_bond_count=row["bridge_bond_count"],
                closure=f"{closure_rmsd:.3f}" if isinstance(closure_rmsd, (float, int)) else "",
                n_clashes=row["n_clashes"],
                worst=row["worst_clash_ratio"] if row["worst_clash_ratio"] is not None else "",
                clash_sources=row["clash_source_summary"] or "",
                topology_n_edges=row["topology_n_edges"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_trace(manifest_path: Path, out_dir: Path, candidate_ids: list[str],
              extra_lattice_scale: float = 1.0) -> list[dict[str, Any]]:
    manifest = _read_manifest(manifest_path)
    missing = [cid for cid in candidate_ids if cid not in manifest]
    if missing:
        raise KeyError(f"Missing candidate ids in manifest: {missing}")

    if out_dir.exists():
        shutil.rmtree(out_dir)
    cif_dir = out_dir / "cifs"
    cif_dir.mkdir(parents=True, exist_ok=True)

    graph_data = enrich_graph_data_with_topology(load_graph("data/ffc/ffc_graph.json"))
    assembler = StructureAssembler(output_dir=str(cif_dir))
    if abs(extra_lattice_scale - 1.0) > 1e-9:
        base_lattice_scale = assembler._compute_lattice_scale

        def _scaled_lattice_scale(graph_data_arg, structure_map_arg, forced_edge_arg=None):
            return base_lattice_scale(
                graph_data_arg, structure_map_arg, forced_edge_arg
            ) * extra_lattice_scale

        assembler._compute_lattice_scale = _scaled_lattice_scale
    rows = []
    for cid in candidate_ids:
        row = manifest[cid]
        coloring = json.loads(row["coloring_json"])
        structure_map = {str(i): name for i, name in enumerate(coloring)}
        cif_path = cif_dir / f"{cid}.cif"
        metadata = assembler.assemble_and_save(
            str(cif_path),
            structure_map,
            graph_data,
            forced_edge=FORCED_EDGE,
            return_metadata=True,
        )
        rows.append(_trace_row(row, metadata, cif_path))

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "regression_ids.json").write_text(
        json.dumps(candidate_ids, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (out_dir / "trace.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    _write_csv(out_dir / "trace.csv", rows)
    _write_markdown(out_dir / "trace.md", rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--candidate-id", action="append", dest="candidate_ids",
                        help="Candidate id to trace. Can be passed multiple times.")
    parser.add_argument("--extra-lattice-scale", type=float, default=1.0,
                        help="Diagnostic-only multiplier applied after assembler lattice scaling.")
    args = parser.parse_args()

    candidate_ids = args.candidate_ids or REGRESSION_IDS
    rows = run_trace(args.manifest, args.out_dir, candidate_ids, args.extra_lattice_scale)
    print(f"Wrote {len(rows)} regression traces to {args.out_dir}")


if __name__ == "__main__":
    main()
