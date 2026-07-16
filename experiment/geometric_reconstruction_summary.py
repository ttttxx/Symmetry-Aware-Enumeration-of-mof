import argparse
import csv
import json
import os
import random
import sys
from collections import Counter, defaultdict
from statistics import median

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from assembler import StructureAssembler
from dft_validation.sanitize import sanitize_one
from main import (
    calculate_topology_requirements,
    enrich_graph_data_with_topology,
    load_graph,
)
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator
from structure_validator import StructureValidator

RANDOM_SEED = 20260703
DEFAULT_PER_CLASS = 5
DEFAULT_RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "reconstruction_rmsd_2026-07-08")

CONFIGS = [
    {
        "topology": "ffc",
        "user_selected_ligands": {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_B_Ch", "4c_Cd_1_Ch"]},
        "forced_edge": "1B_4H_Ch",
        "sampling": {"mode": "stratified", "per_class": DEFAULT_PER_CLASS},
    },
    {
        "topology": "fjh",
        "user_selected_ligands": {3: ["3c_B_Ch"], 4: ["4c_B_Ch", "4c_Cd_1_Ch"]},
        "forced_edge": "1B_4H_Ch",
        "sampling": {"mode": "stratified", "per_class": DEFAULT_PER_CLASS},
    },
    {
        "topology": "vba",
        "user_selected_ligands": {4: ["4c_B_Ch", "4c_Cd_1_Ch", "4c_Cu_1_Ch"]},
        "forced_edge": "1B_4H_Ch",
        "sampling": {"mode": "stratified", "per_class": DEFAULT_PER_CLASS},
    },
    {
        "topology": "hwx",
        "user_selected_ligands": {4: ["4c_Cu_1_Ch"], 6: ["6c_Cu_1_Ch"]},
        "forced_edge": "BDC",
        "sampling": {"mode": "stratified", "per_class": DEFAULT_PER_CLASS},
    },
    {
        "topology": "nbo",
        "user_selected_ligands": {4: ["4c_Cu_1_Ch"]},
        "forced_edge": "BDC",
        "sampling": {"mode": "stratified", "per_class": DEFAULT_PER_CLASS},
    },
]

def _graph_path(topology: str) -> str:
    return os.path.join(PROJECT_ROOT, "data", topology, f"{topology}_graph.json")

def _make_output_dir(base_dir: str, topology: str) -> str:
    path = os.path.join(base_dir, topology)
    os.makedirs(path, exist_ok=True)
    return path

def _classify_sym_classes(results):
    all_stabs = sorted(set(r["stabilizer_size"] for r in results), reverse=True)
    high_vals, med_vals = set(), set()
    if len(all_stabs) <= 2:
        if all_stabs and all_stabs[0] > 1:
            high_vals.add(all_stabs[0])
    elif len(all_stabs) == 3:
        high_vals.add(all_stabs[0])
        if all_stabs[1] > 1:
            med_vals.add(all_stabs[1])
    else:
        high_vals.update(all_stabs[:2])
        for v in all_stabs[2:]:
            if v > 1:
                med_vals.add(v)

    classified = []
    for r in results:
        s = r["stabilizer_size"]
        sym_class = "High" if s in high_vals else ("Med" if s in med_vals else "Low")
        classified.append({**r, "sym_class": sym_class})
    return classified

def _enumerate_case(config):
    topology = config["topology"]
    graph_data = enrich_graph_data_with_topology(load_graph(_graph_path(topology)))
    node_degrees = [n.get("degree", 0) for n in graph_data["nodes"]]
    topology_requirements = calculate_topology_requirements(graph_data)
    perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
    all_ligands = sorted(
        set(lig for ligs in config["user_selected_ligands"].values() for lig in ligs)
    )
    results = MOFEnumerator(
        perms,
        graph_data["n_nodes"],
        verbose=False,
        topology_requirements=topology_requirements,
        user_selected_ligands=config["user_selected_ligands"],
        node_degrees=node_degrees,
    ).enumerate(all_ligands)
    return graph_data, perms, _classify_sym_classes(results)

def _sample_case(classified_results, gamma_size, sampling, seed):
    rng = random.Random(seed)
    mode = sampling.get("mode", "stratified")
    if mode != "stratified":
        raise ValueError(f"Unsupported sampling mode: {mode}")

    per_class = int(sampling.get("per_class", DEFAULT_PER_CLASS))
    grouped = defaultdict(list)
    for r in classified_results:
        grouped[r["sym_class"]].append(r)

    sampled = []
    for sym_class in ("High", "Med", "Low"):
        pool = sorted(grouped.get(sym_class, []), key=lambda x: x["id"])
        if len(pool) <= per_class:
            chosen = pool
        else:
            chosen = rng.sample(pool, per_class)
            chosen = sorted(chosen, key=lambda x: x["id"])
        for r in chosen:
            sampled.append(
                {
                    "result": r,
                    "sym_class": sym_class,
                    "rho": r["stabilizer_size"] / gamma_size,
                }
            )
    return sampled

def _safe_stat_median(values):
    vals = [v for v in values if v is not None]
    return median(vals) if vals else None

def _safe_stat_mean(values):
    vals = [v for v in values if v is not None]
    return sum(vals) / len(vals) if vals else None

def _safe_stat_percentile(values, percentile):
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return None
    if len(vals) == 1:
        return vals[0]
    rank = (len(vals) - 1) * percentile / 100.0
    lo = int(rank)
    hi = min(lo + 1, len(vals) - 1)
    frac = rank - lo
    return vals[lo] * (1.0 - frac) + vals[hi] * frac

def _safe_closure_rmsd(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return None
    return (sum(v * v for v in vals) / len(vals)) ** 0.5

def _format_float(value, digits=3):
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"

def _count_str(numerator, denominator):
    if denominator is None:
        return f"{numerator}/NA"
    return f"{numerator}/{denominator}"

def _format_percent(value, digits=1):
    if value is None:
        return "NA"
    return f"{value * 100:.{digits}f}\\%"

def _format_error_list(values):
    if not values:
        return ""
    return ";".join(f"{float(v):.6f}" for v in values)

def _parse_error_list(value):
    if value in (None, ""):
        return []
    if isinstance(value, (list, tuple)):
        return [float(v) for v in value]
    return [float(v) for v in str(value).split(";") if v]

def _closure_summary_from_rows(rows):
    errors = []
    metric_valid = 0
    for row in rows:
        row_errors = _parse_error_list(row.get("closure_errors_A"))
        if row_errors:
            metric_valid += 1
            errors.extend(row_errors)

    return {
        "metric_valid": metric_valid,
        "n_closure_pairs": len(errors),
        "closure_rmsd_A": _safe_closure_rmsd(errors),
        "closure_median_A": _safe_stat_median(errors),
        "closure_p95_A": _safe_stat_percentile(errors, 95),
        "closure_max_A": max(errors) if errors else None,
        "within_0.8A_fraction": (
            sum(1 for d in errors if d <= 0.8) / len(errors)
            if errors else None
        ),
    }

def _assemble_and_audit_one(assembler, graph_data, forced_edge, base_dir, topology, rec, idx):
    result = rec["result"]
    structure_map = {str(i): color for i, color in enumerate(result["coloring"])}
    cif_name = f"{topology}_{rec['sym_class']}_rho{rec['rho']:.4f}_{idx:03d}.cif"
    cif_path = os.path.join(base_dir, cif_name)

    row = {
        "topology": topology,
        "structure_id": result["id"],
        "sym_class": rec["sym_class"],
        "rho": rec["rho"],
        "stabilizer_size": result["stabilizer_size"],
        "assembly_ok": False,
        "assembly_error": "",
        "sanitize_status": "",
        "sanitize_ok": False,
        "n_atoms": None,
        "scale_factor": None,
        "merge_count": None,
        "bridge_bond_count": None,
        "closure_metric_ok": False,
        "n_closure_pairs": 0,
        "closure_rmsd_A": None,
        "closure_median_A": None,
        "closure_p95_A": None,
        "closure_max_A": None,
        "within_0.8A_fraction": None,
        "closure_errors_A": "",
        "validation_passed": None,
        "clash_pass": None,
        "clash_n": None,
        "clash_min_ratio": None,
        "topology_pass": None,
        "topology_n_assembled": None,
        "topology_n_target": None,
        "topology_e_assembled": None,
        "topology_e_target": None,
        "symmetry_pass": None,
        "detected_sg": "",
        "symmetry_quality": "",
        "symmetry_n_ops": None,
        "output_cif": cif_path,
    }

    try:
        metadata = assembler.assemble_and_save(
            filename=cif_path,
            structure_map=structure_map,
            graph_data=graph_data,
            forced_edge=forced_edge,
            return_metadata=True,
        )
        row["assembly_ok"] = True
        row["n_atoms"] = metadata.get("n_atoms")
        row["scale_factor"] = metadata.get("scale_factor")
        row["merge_count"] = metadata.get("merge_count")
        row["bridge_bond_count"] = metadata.get("bridge_bond_count")
        row["validation_passed"] = metadata.get("validation_passed")
        closure = metadata.get("connection_closure", {})
        row["n_closure_pairs"] = closure.get("n_closure_pairs", 0)
        row["closure_rmsd_A"] = closure.get("closure_rmsd_A")
        row["closure_median_A"] = closure.get("closure_median_A")
        row["closure_p95_A"] = closure.get("closure_p95_A")
        row["closure_max_A"] = closure.get("closure_max_A")
        row["within_0.8A_fraction"] = closure.get("within_0.8A_fraction")
        row["closure_errors_A"] = _format_error_list(closure.get("closure_errors_A", []))
        row["closure_metric_ok"] = row["n_closure_pairs"] > 0 and row["closure_rmsd_A"] is not None

        validation = metadata.get("validation", {})
        clash = validation.get("clash", {})
        topo = validation.get("topology", {})
        sym = validation.get("symmetry", {})

        row["clash_pass"] = clash.get("passed")
        row["clash_n"] = clash.get("n_clashes")
        row["clash_min_ratio"] = clash.get("min_ratio")

        row["topology_pass"] = topo.get("passed")
        row["topology_n_assembled"] = topo.get("n_assembled")
        row["topology_n_target"] = topo.get("n_tg")
        row["topology_e_assembled"] = topo.get("e_assembled")
        row["topology_e_target"] = topo.get("e_tg")

        row["symmetry_pass"] = sym.get("passed")
        row["detected_sg"] = sym.get("detected_sg", "")
        row["symmetry_quality"] = sym.get("symmetry_quality", "")
        row["symmetry_n_ops"] = sym.get("n_symops")

        sanitize_record = {
            "id": f"{topology}_{result['id']}",
            "topology": topology,
            "status": "assembled",
            "init_cif": cif_path,
            "n_atoms": row["n_atoms"],
        }
        sanitized = sanitize_one(sanitize_record, os.path.join(base_dir, "_sanitize"))
        row["sanitize_status"] = sanitized.get("status", "")
        row["sanitize_ok"] = sanitized.get("status") in ("clean", "preclean_done")
    except Exception as e:
        row["assembly_error"] = str(e)

    return row

def _write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

def _summarize_by_topology(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["topology"]].append(row)

    summary_rows = []
    for topology, items in sorted(grouped.items()):
        attempted = len(items)
        assembly_ok = sum(1 for x in items if x["assembly_ok"])
        sanitize_ok = sum(1 for x in items if x["sanitize_ok"])
        topology_pass = sum(1 for x in items if x["topology_pass"] is True)
        clash_pass = sum(1 for x in items if x["clash_pass"] is True)
        relaxed_validation = sum(
            1 for x in items if x["assembly_ok"] and x["topology_pass"] is True and x["clash_pass"] is True
        )
        closure_summary = _closure_summary_from_rows(items)
        summary_rows.append(
            {
                "topology": topology,
                "attempted": attempted,
                "assembly_ok": assembly_ok,
                "metric_valid": closure_summary["metric_valid"],
                "n_closure_pairs": closure_summary["n_closure_pairs"],
                "closure_rmsd_A": closure_summary["closure_rmsd_A"],
                "closure_median_A": closure_summary["closure_median_A"],
                "closure_p95_A": closure_summary["closure_p95_A"],
                "closure_max_A": closure_summary["closure_max_A"],
                "within_0.8A_fraction": closure_summary["within_0.8A_fraction"],
                "assembly_success_rate": assembly_ok / attempted if attempted else None,
                "sanitize_ok": sanitize_ok,
                "sanitize_success_rate": sanitize_ok / assembly_ok if assembly_ok else None,
                "topology_pass": topology_pass,
                "topology_recovery_rate": topology_pass / assembly_ok if assembly_ok else None,
                "clash_pass": clash_pass,
                "clash_free_rate": clash_pass / assembly_ok if assembly_ok else None,
                "relaxed_validation_pass": relaxed_validation,
                "relaxed_validation_rate": relaxed_validation / assembly_ok if assembly_ok else None,
                "median_clash_min_ratio": _safe_stat_median([x["clash_min_ratio"] for x in items if x["assembly_ok"]]),
                "mean_clash_count": _safe_stat_mean([x["clash_n"] for x in items if x["assembly_ok"]]),
            }
        )
    return summary_rows

def _summarize_overall(rows, by_topology):
    attempted = len(rows)
    assembly_ok = sum(1 for x in rows if x["assembly_ok"])
    sanitize_ok = sum(1 for x in rows if x["sanitize_ok"])
    topology_pass = sum(1 for x in rows if x["topology_pass"] is True)
    clash_pass = sum(1 for x in rows if x["clash_pass"] is True)
    relaxed_validation = sum(
        1 for x in rows if x["assembly_ok"] and x["topology_pass"] is True and x["clash_pass"] is True
    )
    sg_counter = Counter(x["detected_sg"] for x in rows if x["assembly_ok"] and x["detected_sg"])
    closure_summary = _closure_summary_from_rows(rows)

    return {
        "attempted": attempted,
        "assembly_ok": assembly_ok,
        "metric_valid": closure_summary["metric_valid"],
        "n_closure_pairs": closure_summary["n_closure_pairs"],
        "closure_rmsd_A": closure_summary["closure_rmsd_A"],
        "closure_median_A": closure_summary["closure_median_A"],
        "closure_p95_A": closure_summary["closure_p95_A"],
        "closure_max_A": closure_summary["closure_max_A"],
        "within_0.8A_fraction": closure_summary["within_0.8A_fraction"],
        "assembly_success_rate": assembly_ok / attempted if attempted else None,
        "sanitize_ok": sanitize_ok,
        "sanitize_success_rate": sanitize_ok / assembly_ok if assembly_ok else None,
        "topology_pass": topology_pass,
        "topology_recovery_rate": topology_pass / assembly_ok if assembly_ok else None,
        "clash_pass": clash_pass,
        "clash_free_rate": clash_pass / assembly_ok if assembly_ok else None,
        "relaxed_validation_pass": relaxed_validation,
        "relaxed_validation_rate": relaxed_validation / assembly_ok if assembly_ok else None,
        "median_clash_min_ratio": _safe_stat_median([x["clash_min_ratio"] for x in rows if x["assembly_ok"]]),
        "mean_clash_count": _safe_stat_mean([x["clash_n"] for x in rows if x["assembly_ok"]]),
        "spacegroup_distribution": dict(sorted(sg_counter.items())),
        "topologies": [row["topology"] for row in by_topology],
    }

def _write_markdown_report(path, overall, by_topology):
    lines = [
        "# Geometric Reconstruction Summary",
        "",
        "## Overall",
        "",
        f"- Attempted reconstructions: {overall['attempted']}",
        f"- Assembly success: {overall['assembly_ok']}/{overall['attempted']} ({_format_float(overall['assembly_success_rate'] * 100 if overall['assembly_success_rate'] is not None else None, 1)}%)",
        f"- Metric-valid reconstructions: {overall['metric_valid']}/{overall['assembly_ok']}",
        f"- Connection-site closure RMSD: {_format_float(overall['closure_rmsd_A'])} Å",
        f"- Connection-site closure median / p95 / max: {_format_float(overall['closure_median_A'])} / {_format_float(overall['closure_p95_A'])} / {_format_float(overall['closure_max_A'])} Å",
        f"- Closure pairs within 0.8 Å: {_format_float(overall['within_0.8A_fraction'] * 100 if overall['within_0.8A_fraction'] is not None else None, 1)}%",
        f"- Sanitizable CIFs: {overall['sanitize_ok']}/{overall['assembly_ok']} ({_format_float(overall['sanitize_success_rate'] * 100 if overall['sanitize_success_rate'] is not None else None, 1)}%)",
        f"- Topology-consistent reconstructions: {overall['topology_pass']}/{overall['assembly_ok']} ({_format_float(overall['topology_recovery_rate'] * 100 if overall['topology_recovery_rate'] is not None else None, 1)}%)",
        f"- Clash-free reconstructions: {overall['clash_pass']}/{overall['assembly_ok']} ({_format_float(overall['clash_free_rate'] * 100 if overall['clash_free_rate'] is not None else None, 1)}%)",
        f"- Combined topology-pass and clash-pass: {overall['relaxed_validation_pass']}/{overall['assembly_ok']} ({_format_float(overall['relaxed_validation_rate'] * 100 if overall['relaxed_validation_rate'] is not None else None, 1)}%)",
        f"- Median minimum vdW ratio: {_format_float(overall['median_clash_min_ratio'])}",
        f"- Mean clash count: {_format_float(overall['mean_clash_count'])}",
        f"- Space-group distribution among assembled structures: `{overall['spacegroup_distribution']}`",
        "",
        "## By Topology",
        "",
        "| Topology | Sampled | Metric-valid | Closure pairs | RMSD (Å) | Median (Å) | p95 (Å) | Max (Å) | ≤0.8 Å |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in by_topology:
        lines.append(
            f"| {row['topology']} | {row['attempted']} | "
            f"{_count_str(row['metric_valid'], row['assembly_ok'])} | "
            f"{row['n_closure_pairs']} | "
            f"{_format_float(row['closure_rmsd_A'])} | "
            f"{_format_float(row['closure_median_A'])} | "
            f"{_format_float(row['closure_p95_A'])} | "
            f"{_format_float(row['closure_max_A'])} | "
            f"{_format_float(row['within_0.8A_fraction'] * 100 if row['within_0.8A_fraction'] is not None else None, 1)}% |"
        )
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

def _write_latex_table(path, by_topology):
    lines = [
        "\\begin{tabular}{lcccccccc}",
        "\\toprule",
        "Topology & Sampled & Metric-valid & Pairs & RMSD (\\AA) & Median (\\AA) & p95 (\\AA) & Max (\\AA) & $\\le$0.8\\,\\AA \\\\",
        "\\midrule",
    ]
    for row in by_topology:
        lines.append(
            f"{row['topology']} & {row['attempted']} & "
            f"{_count_str(row['metric_valid'], row['assembly_ok'])} & "
            f"{row['n_closure_pairs']} & "
            f"{_format_float(row['closure_rmsd_A'])} & "
            f"{_format_float(row['closure_median_A'])} & "
            f"{_format_float(row['closure_p95_A'])} & "
            f"{_format_float(row['closure_max_A'])} & "
            f"{_format_percent(row['within_0.8A_fraction'])} \\\\"
        )
    lines.extend(["\\bottomrule", "\\end{tabular}"])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

def run_summary(configs=None, base_dir=DEFAULT_RESULTS_DIR, seed=RANDOM_SEED):
    configs = configs or CONFIGS
    os.makedirs(base_dir, exist_ok=True)
    all_rows = []

    for config_idx, config in enumerate(configs):
        topology = config["topology"]
        print(f"\n=== Reconstruction summary: {topology} ===")
        graph_data, perms, classified = _enumerate_case(config)
        gamma_size = len(perms)
        sampled = _sample_case(
            classified,
            gamma_size=gamma_size,
            sampling=config.get("sampling", {}),
            seed=seed + config_idx,
        )
        print(f"  Enumerated {len(classified)} symmetry-distinct structures; sampled {len(sampled)} for audit.")
        output_dir = _make_output_dir(base_dir, topology)
        assembler = StructureAssembler(
            output_dir=output_dir,
            validator=StructureValidator(strict_mode=False, check_symmetry=False),
            run_symmetry=False,
        )

        for idx, rec in enumerate(sampled):
            row = _assemble_and_audit_one(
                assembler=assembler,
                graph_data=graph_data,
                forced_edge=config["forced_edge"],
                base_dir=output_dir,
                topology=topology,
                rec=rec,
                idx=idx,
            )
            all_rows.append(row)

    records_csv = os.path.join(base_dir, "reconstruction_records.csv")
    _write_csv(records_csv, all_rows)

    by_topology = _summarize_by_topology(all_rows)
    topology_csv = os.path.join(base_dir, "reconstruction_by_topology.csv")
    _write_csv(topology_csv, by_topology)

    overall = _summarize_overall(all_rows, by_topology)
    overall_csv = os.path.join(base_dir, "reconstruction_overall.csv")
    _write_csv(overall_csv, [overall])

    report_md = os.path.join(base_dir, "reconstruction_report.md")
    _write_markdown_report(report_md, overall, by_topology)

    table_tex = os.path.join(base_dir, "reconstruction_table_main.tex")
    _write_latex_table(table_tex, by_topology)

    return {
        "records_csv": records_csv,
        "topology_csv": topology_csv,
        "overall_csv": overall_csv,
        "report_md": report_md,
        "table_tex": table_tex,
        "n_records": len(all_rows),
    }

def _select_configs(topologies=None, per_class=None):
    configs = []
    wanted = set(topologies) if topologies else None
    for cfg in CONFIGS:
        if wanted and cfg["topology"] not in wanted:
            continue
        copied = {
            "topology": cfg["topology"],
            "user_selected_ligands": cfg["user_selected_ligands"],
            "forced_edge": cfg["forced_edge"],
            "sampling": dict(cfg.get("sampling", {})),
        }
        if per_class is not None:
            copied["sampling"]["per_class"] = per_class
        configs.append(copied)
    if not configs:
        raise ValueError("No matching topologies found in CONFIGS.")
    return configs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Audit geometric reconstruction quality across representative MOF topologies.")
    parser.add_argument("--topology", action="append", dest="topologies",
                        help="Topology to audit. May be passed multiple times. Defaults to all preset CONFIGS.")
    parser.add_argument("--per-class", type=int, default=None,
                        help="Override the number of sampled structures per symmetry class.")
    parser.add_argument("--base-dir", default=DEFAULT_RESULTS_DIR,
                        help="Output directory for summary artifacts.")
    parser.add_argument("--seed", type=int, default=RANDOM_SEED,
                        help="Random seed for stratified sampling.")
    args = parser.parse_args()

    summary = run_summary(
        configs=_select_configs(topologies=args.topologies, per_class=args.per_class),
        base_dir=args.base_dir,
        seed=args.seed,
    )
    print("\nGenerated outputs:")
    for key, value in summary.items():
        print(f"  {key}: {value}")
