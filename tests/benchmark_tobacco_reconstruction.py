#!/usr/bin/env python
"""Run a non-interactive ToBaCCo reconstruction benchmark.

The benchmark turns the reference CIFs under ``tests/tobacco_output_cifs``
into a reproducible engineering diagnostic:

* parse each ToBaCCo-style case name;
* run the local enumeration/reconstruction pipeline;
* compare generated CIFs against the reference CIF;
* write per-case, per-topology, and Markdown summary reports.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import os
import re
import sys
import time
from collections import Counter, defaultdict
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

from pymatgen.core import Structure

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
os.chdir(PROJECT_ROOT)

from external_lib.node import Node
from main import calculate_topology_requirements, enrich_graph_data_with_topology, load_graph
from tests.comparator import StructureReport, compare_structures
from tests.name_parser import parse_filename


REF_DIR = PROJECT_ROOT / "tests" / "tobacco_output_cifs"
DEFAULT_OUT_DIR = PROJECT_ROOT / "results" / "tobacco_reconstruction_benchmark"
_NODE_X_COUNT_CACHE: dict[str, int | None] = {}


@dataclass
class CaseInput:
    ref_file: str
    topology: str
    nodes: list[str]
    edge: str


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _degree_from_node(node_name: str) -> int | None:
    match = re.match(r"(\d+)c_", node_name)
    if not match:
        return None
    return int(match.group(1))


def _json_dump(value) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _node_x_count(node_name: str) -> int | None:
    if node_name not in _NODE_X_COUNT_CACHE:
        try:
            node = Node(str(PROJECT_ROOT / "nodes" / f"{node_name}.cif"))
            _NODE_X_COUNT_CACHE[node_name] = len(node.calculate_X_cart_vectors(False))
        except Exception:
            _NODE_X_COUNT_CACHE[node_name] = None
    return _NODE_X_COUNT_CACHE[node_name]


def _load_topology_context(topology: str) -> tuple[dict | None, dict[int, int], list[int]]:
    graph_path = PROJECT_ROOT / "data" / topology / f"{topology}_graph.json"
    if not graph_path.exists():
        return None, {}, []
    with redirect_stdout(io.StringIO()):
        graph_data = enrich_graph_data_with_topology(load_graph(str(graph_path)))
    requirements = calculate_topology_requirements(graph_data)
    node_degrees = [node.get("degree", 0) for node in graph_data["nodes"]]
    return graph_data, requirements, node_degrees


def _scan_reference_files(limit: int | None = None, cases: Iterable[str] | None = None) -> list[str]:
    selected = {Path(c).name for c in cases or []}
    files = sorted(p.name for p in REF_DIR.glob("*.cif"))
    if selected:
        files = [f for f in files if f in selected or Path(f).stem in selected]
    if limit is not None:
        files = files[:limit]
    return files


def _parse_case(ref_file: str) -> tuple[CaseInput | None, str | None]:
    try:
        topology, nodes, edge = parse_filename(ref_file)
    except Exception as exc:
        return None, f"parse_error: {exc}"
    return CaseInput(ref_file=ref_file, topology=topology, nodes=nodes, edge=edge), None


def _check_case_inputs(case: CaseInput) -> list[str]:
    issues = []
    graph_path = PROJECT_ROOT / "data" / case.topology / f"{case.topology}_graph.json"
    if not graph_path.exists():
        issues.append("missing_topology_graph")

    edge_path = PROJECT_ROOT / "edges" / f"{case.edge}.cif"
    if not edge_path.exists():
        issues.append("missing_edge_cif")

    for node in case.nodes:
        if not (PROJECT_ROOT / "nodes" / f"{node}.cif").exists():
            issues.append(f"missing_node_cif:{node}")
        if _degree_from_node(node) is None:
            issues.append(f"invalid_node_degree:{node}")

    return issues


def _build_ligand_selection(nodes: list[str]) -> dict[int, list[str]]:
    selected: dict[int, list[str]] = defaultdict(list)
    for node in nodes:
        degree = _degree_from_node(node)
        if degree is None:
            raise ValueError(f"Cannot parse connectivity degree from node name: {node}")
        if node not in selected[degree]:
            selected[degree].append(node)
    return dict(selected)


def _ligand_degree_map(nodes: list[str]) -> dict[int, list[str]]:
    return _build_ligand_selection(nodes)


def _ligand_x_count_map(nodes: list[str]) -> dict[str, int | None]:
    return {node: _node_x_count(node) for node in sorted(set(nodes))}


def _degree_compatibility(requirements: dict[int, int], ligand_map: dict[int, list[str]]) -> tuple[list[int], list[int]]:
    missing = sorted(deg for deg in requirements if deg not in ligand_map or not ligand_map[deg])
    extra = sorted(deg for deg in ligand_map if deg not in requirements)
    return missing, extra


def _status_rank(status: str) -> int:
    if status == StructureReport.STATUS_BAD:
        return 2
    if status == StructureReport.STATUS_WARN:
        return 1
    return 0


def _report_score(report: StructureReport) -> tuple[int, int, int]:
    return (
        len(report.errors),
        len(report.warnings),
        sum(_status_rank(status) for _, status, _ in report.checks),
    )


def _check_detail(report: StructureReport, label: str) -> tuple[str, str]:
    for check_label, status, detail in report.checks:
        if check_label == label:
            return status, detail
    return "not_checked", ""


def _structure_metrics(ref_path: Path, gen_path: Path | None) -> dict[str, str]:
    metrics = {
        "ref_atoms": "",
        "gen_atoms": "",
        "atom_ratio": "",
        "ref_volume": "",
        "gen_volume": "",
        "volume_ratio": "",
    }
    try:
        ref = Structure.from_file(ref_path)
        metrics["ref_atoms"] = str(len(ref))
        metrics["ref_volume"] = f"{ref.lattice.volume:.3f}"
    except Exception:
        return metrics

    if gen_path is None:
        return metrics

    try:
        gen = Structure.from_file(gen_path)
    except Exception:
        return metrics

    metrics["gen_atoms"] = str(len(gen))
    if len(ref):
        metrics["atom_ratio"] = f"{len(gen) / len(ref):.4f}"
    metrics["gen_volume"] = f"{gen.lattice.volume:.3f}"
    if ref.lattice.volume > 0:
        metrics["volume_ratio"] = f"{gen.lattice.volume / ref.lattice.volume:.4f}"
    return metrics


def _compare_generated(ref_path: Path, cif_paths: list[str]) -> tuple[Path | None, StructureReport | None]:
    if not cif_paths:
        return None, None
    scored = []
    for cif_path in cif_paths:
        report = compare_structures(str(ref_path), cif_path)
        scored.append((_report_score(report), Path(cif_path), report))
    scored.sort(key=lambda item: item[0])
    _, best_path, best_report = scored[0]
    return best_path, best_report


def _failure_class(input_issues: list[str], error: str, report: StructureReport | None,
                   diagnostic: str = "", degree_leaks: str = "") -> str:
    if input_issues:
        if "missing_topology_graph" in input_issues:
            return "input_missing_topology"
        if "input_degree_incompatible" in input_issues:
            return "input_degree_incompatible"
        if any(i.startswith("missing_node_cif") for i in input_issues):
            return "input_missing_node"
        if "missing_edge_cif" in input_issues:
            return "input_missing_edge"
        return "input_invalid"

    if error:
        if diagnostic.startswith("shape_mismatch"):
            if degree_leaks:
                return "enumerator_degree_leak"
            return "assembler_anchor_mismatch"
        if "degree_incompatible" in error:
            return "input_degree_incompatible"
        if "enumeration_no_structures" in error:
            return "enumeration_no_structures"
        if diagnostic.startswith("shape_mismatch"):
            return "assembly_shape_mismatch"
        if "no_cif_generated" in error:
            return "assembly_no_cif"
        return "assembly_exception"

    if report is None:
        return "not_evaluated"

    if report.errors:
        return "comparison_severe_difference"
    if report.warnings:
        return "comparison_warning_difference"
    return "passed_all_checks"


def _extract_log_diagnostic(log_path: Path) -> str:
    if not log_path.exists():
        return ""
    text = log_path.read_text(errors="replace")
    match = re.search(r"配体连接点数 \((\d+)\) 与拓扑锚点数 \((\d+)\) 不匹配", text)
    if match:
        return f"shape_mismatch: ligand_connectors={match.group(1)}, topology_anchors={match.group(2)}"
    if "[Validation] WARNING" in text:
        failed = re.findall(r"\[FAIL\] ([^:]+): ([^\n]+)", text)
        if failed:
            return "; ".join(f"validation_{name.strip()}: {detail.strip()}" for name, detail in failed)
        return "validation_warning"
    if "[Error]" in text:
        err_lines = [line.strip() for line in text.splitlines() if "[Error]" in line]
        return err_lines[-1] if err_lines else "assembly_error"
    return ""


def _extract_prealign_trace(log_path: Path) -> list[dict[str, str]]:
    if not log_path.exists():
        return []
    trace = []
    pattern = re.compile(
        r"\[PreAlign\]\s+node=(?P<node>\S+)\s+cif=(?P<cif>\S+)\s+"
        r"node_x_count=(?P<node_x_count>\S+)\s+topology_anchor_count=(?P<topology_anchor_count>\S+)"
    )
    for line in log_path.read_text(errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            trace.append(match.groupdict())
    return trace


def _extract_atom_trace(log_path: Path) -> dict[str, str]:
    if not log_path.exists():
        return {}
    for line in log_path.read_text(errors="replace").splitlines():
        if "[AtomTrace]" not in line:
            continue
        values = {}
        for part in line.split("[AtomTrace]", 1)[1].strip().split():
            if "=" in part:
                key, value = part.split("=", 1)
                values[key] = value
        return values
    return {}


def _coloring_degree_trace(graph_data: dict | None, coloring: list[str]) -> list[dict[str, object]]:
    if graph_data is None or not coloring:
        return []
    trace = []
    for idx, ligand in enumerate(coloring):
        node = graph_data["nodes"][idx]
        graph_degree = node.get("degree", 0)
        ligand_x = _node_x_count(ligand)
        trace.append({
            "node_index": idx,
            "node_id": str(node.get("id", idx)),
            "graph_degree": graph_degree,
            "ligand": ligand,
            "ligand_declared_degree": _degree_from_node(ligand),
            "ligand_x_count": ligand_x,
            "degree_match": ligand_x == graph_degree,
        })
    return trace


def _degree_leaks_from_trace(trace: list[dict[str, object]]) -> list[dict[str, object]]:
    return [entry for entry in trace if entry.get("degree_match") is False]


def _run_one_case(case: CaseInput, out_dir: Path, max_cifs: int, max_results: int) -> dict[str, str]:
    ref_path = REF_DIR / case.ref_file
    case_name = Path(case.ref_file).stem
    log_path = out_dir / "logs" / f"{_safe_name(case_name)}.log"

    row = {
        "case": case_name,
        "ref_file": case.ref_file,
        "topology": case.topology,
        "nodes": ";".join(case.nodes),
        "edge": case.edge,
        "input_issues": "",
        "topology_requirements": "",
        "ligand_degree_map": "",
        "ligand_x_counts": "",
        "missing_topology_degrees": "",
        "extra_ligand_degrees": "",
        "first_coloring": "",
        "coloring_degree_trace": "",
        "coloring_degree_leaks": "",
        "prealign_trace": "",
        "atom_trace": "",
        "generated_cif": "",
        "n_generated_cifs": "0",
        "failure_class": "",
        "errors": "",
        "warnings": "",
        "error_message": "",
        "diagnostic": "",
        "log_file": str(log_path.relative_to(out_dir)),
    }

    input_issues = _check_case_inputs(case)
    row["input_issues"] = ";".join(input_issues)
    if input_issues:
        row["failure_class"] = _failure_class(input_issues, "", None)
        row.update(_structure_metrics(ref_path, None))
        return row

    try:
        user_selected_ligands = _build_ligand_selection(case.nodes)
    except Exception as exc:
        row["error_message"] = str(exc)
        row["failure_class"] = "input_invalid"
        row.update(_structure_metrics(ref_path, None))
        return row

    graph_data, requirements, _node_degrees = _load_topology_context(case.topology)
    ligand_map = _ligand_degree_map(case.nodes)
    missing_degrees, extra_degrees = _degree_compatibility(requirements, ligand_map)
    row["topology_requirements"] = _json_dump(requirements)
    row["ligand_degree_map"] = _json_dump(ligand_map)
    row["ligand_x_counts"] = _json_dump(_ligand_x_count_map(case.nodes))
    row["missing_topology_degrees"] = ";".join(str(d) for d in missing_degrees)
    row["extra_ligand_degrees"] = ";".join(str(d) for d in extra_degrees)

    if missing_degrees:
        row["input_issues"] = ";".join(filter(None, [row["input_issues"], "input_degree_incompatible"]))
        row["error_message"] = (
            "degree_incompatible: missing topology degrees "
            f"{missing_degrees}; ligand degrees {sorted(ligand_map)}"
        )
        row["failure_class"] = "input_degree_incompatible"
        row.update(_structure_metrics(ref_path, None))
        return row

    cif_paths: list[str] = []
    final_output = {}
    error_message = ""
    root_logger = logging.getLogger()
    old_handlers = list(root_logger.handlers)
    old_level = root_logger.level

    with open(log_path, "w") as log_fh:
        for handler in old_handlers:
            root_logger.removeHandler(handler)
        log_handler = logging.StreamHandler(log_fh)
        log_handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s:%(message)s"))
        root_logger.addHandler(log_handler)
        root_logger.setLevel(logging.WARNING)

        print(f"[Benchmark] {case.ref_file}", file=log_fh)
        print(f"[Benchmark] topology={case.topology}", file=log_fh)
        print(f"[Benchmark] nodes={case.nodes}", file=log_fh)
        print(f"[Benchmark] edge={case.edge}", file=log_fh)
        try:
            from main import run_enumeration

            with redirect_stdout(log_fh), redirect_stderr(log_fh):
                cif_paths, _, final_output = run_enumeration(
                    data_prefix=case.topology,
                    user_selected_ligands=user_selected_ligands,
                    forced_edge=case.edge,
                    max_cifs=max_cifs,
                    max_results=max_results,
                )
        except Exception as exc:
            error_message = str(exc)
            print(f"[Benchmark] exception={error_message}", file=log_fh)
        finally:
            root_logger.removeHandler(log_handler)
            for handler in old_handlers:
                root_logger.addHandler(handler)
            root_logger.setLevel(old_level)

    diagnostic = _extract_log_diagnostic(log_path)
    row["diagnostic"] = diagnostic
    row["prealign_trace"] = _json_dump(_extract_prealign_trace(log_path))
    row["atom_trace"] = _json_dump(_extract_atom_trace(log_path))

    structures = final_output.get("structures", []) if isinstance(final_output, dict) else []
    if structures:
        first_coloring = structures[0].get("coloring", [])
        row["first_coloring"] = _json_dump(first_coloring)
        degree_trace = _coloring_degree_trace(graph_data, first_coloring)
        degree_leaks = _degree_leaks_from_trace(degree_trace)
        row["coloring_degree_trace"] = _json_dump(degree_trace)
        row["coloring_degree_leaks"] = _json_dump(degree_leaks)

    row["n_generated_cifs"] = str(len(cif_paths))
    if error_message:
        row["error_message"] = error_message
        row["failure_class"] = _failure_class([], error_message, None, diagnostic, row["coloring_degree_leaks"])
        row.update(_structure_metrics(ref_path, None))
        return row

    if not structures:
        row["error_message"] = "enumeration_no_structures"
        row["failure_class"] = "enumeration_no_structures"
        row.update(_structure_metrics(ref_path, None))
        return row

    if not cif_paths:
        row["error_message"] = "no_cif_generated"
        row["failure_class"] = _failure_class([], "no_cif_generated", None, diagnostic, row["coloring_degree_leaks"])
        row.update(_structure_metrics(ref_path, None))
        return row

    best_path, report = _compare_generated(ref_path, cif_paths)
    row["generated_cif"] = str(best_path) if best_path else ""
    row["failure_class"] = _failure_class([], "", report, diagnostic, row["coloring_degree_leaks"])
    if report is not None:
        row["errors"] = ";".join(report.errors)
        row["warnings"] = ";".join(report.warnings)
        row["error_details"] = _json_dump({
            label: detail for label, status, detail in report.checks
            if status == StructureReport.STATUS_BAD
        })
        row["warning_details"] = _json_dump({
            label: detail for label, status, detail in report.checks
            if status == StructureReport.STATUS_WARN
        })
        for label in ("加载参考CIF", "加载生成CIF", "总原子数", "晶胞参数", "晶胞体积", "空间群", "元素组成",
                      "参考结构原子间距", "生成结构原子间距", "生成结构成键间距"):
            status, detail = _check_detail(report, label)
            key = {
                "加载参考CIF": "load_ref_check",
                "加载生成CIF": "load_gen_check",
                "总原子数": "atom_count_check",
                "晶胞参数": "lattice_parameter_check",
                "晶胞体积": "volume_check",
                "空间群": "space_group_check",
                "元素组成": "composition_check",
                "参考结构原子间距": "ref_min_distance_check",
                "生成结构原子间距": "gen_min_distance_check",
                "生成结构成键间距": "gen_bond_distance_check",
            }[label]
            row[key] = status
            row[f"{key}_detail"] = detail

    row.update(_structure_metrics(ref_path, best_path))
    return row


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _topology_summary(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["topology"]].append(row)

    summary = []
    for topology, group in sorted(grouped.items()):
        classes = Counter(r["failure_class"] for r in group)
        generated = sum(int(r.get("n_generated_cifs") or 0) > 0 for r in group)
        severe = classes["comparison_severe_difference"]
        warnings = classes["comparison_warning_difference"]
        passed = classes["passed_all_checks"]
        summary.append({
            "topology": topology,
            "cases": str(len(group)),
            "generated_cases": str(generated),
            "passed_all_checks": str(passed),
            "warning_difference": str(warnings),
            "severe_difference": str(severe),
            "input_missing_topology": str(classes["input_missing_topology"]),
            "input_degree_incompatible": str(classes["input_degree_incompatible"]),
            "enumeration_no_structures": str(classes["enumeration_no_structures"]),
            "enumerator_degree_leak": str(classes["enumerator_degree_leak"]),
            "assembler_anchor_mismatch": str(classes["assembler_anchor_mismatch"]),
            "shape_mismatch": str(classes["assembly_shape_mismatch"]),
            "assembly_no_cif": str(classes["assembly_no_cif"]),
            "assembly_exception": str(classes["assembly_exception"]),
            "main_failure_classes": ";".join(
                f"{name}:{count}" for name, count in classes.most_common() if count
            ),
        })
    return summary


def _write_markdown(path: Path, rows: list[dict[str, str]], summary_rows: list[dict[str, str]],
                    elapsed_s: float, args: argparse.Namespace) -> None:
    total = len(rows)
    classes = Counter(r["failure_class"] for r in rows)
    generated = sum(int(r.get("n_generated_cifs") or 0) > 0 for r in rows)
    lines = [
        "# ToBaCCo Reconstruction Benchmark",
        "",
        f"- Generated at: `{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}`",
        f"- Reference directory: `{REF_DIR.relative_to(PROJECT_ROOT)}`",
        f"- Cases attempted: `{total}`",
        f"- Cases with generated CIF: `{generated}`",
        f"- Max CIFs per case: `{args.max_cifs}`",
        f"- Max enumeration results per case: `{args.max_results}`",
        f"- Runtime: `{elapsed_s:.1f} s`",
        "",
        "## Outcome Classes",
        "",
        "| Class | Count |",
        "|---|---:|",
    ]
    for name, count in classes.most_common():
        lines.append(f"| `{name}` | {count} |")

    lines.extend([
        "",
        "## Topology Summary",
        "",
        "| Topology | Cases | Generated | Passed | Warning diff. | Severe diff. | Missing topology | Input degree incompatible | No structures | Degree leak | Anchor mismatch | No CIF | Exception |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for row in summary_rows:
        lines.append(
            f"| `{row['topology']}` | {row['cases']} | {row['generated_cases']} | "
            f"{row['passed_all_checks']} | {row['warning_difference']} | "
            f"{row['severe_difference']} | {row['input_missing_topology']} | "
            f"{row['input_degree_incompatible']} | {row['enumeration_no_structures']} | "
            f"{row['enumerator_degree_leak']} | {row['assembler_anchor_mismatch']} | "
            f"{row['assembly_no_cif']} | {row['assembly_exception']} |"
        )

    lines.extend([
        "",
        "## Engineering Interpretation",
        "",
        "- `input_missing_topology` means the ToBaCCo reference exists but this repository currently lacks the matching topology graph input.",
        "- `input_degree_incompatible` means the parsed ToBaCCo case does not provide ligands for at least one topology node degree, so reconstruction is skipped before enumeration.",
        "- `enumeration_no_structures` means strict degree constraints left no valid coloring.",
        "- `enumerator_degree_leak` means a wrong-degree ligand reached assembler despite strict constraints and should be treated as an enumerator bug.",
        "- `assembler_anchor_mismatch` means coloring degrees look valid but the pre-alignment anchor count still differs from the node CIF X count.",
        "- `assembly_no_cif` means enumeration ran but no reconstructed CIF survived the current rigid-body reconstruction path.",
        "- `comparison_warning_difference` marks generated structures that are loadable but differ from the ToBaCCo reference in cell, composition, space group, or distance diagnostics.",
        "- `comparison_severe_difference` marks generated structures with at least one severe comparator finding.",
        "",
        "## Next Debug Targets",
        "",
    ])

    debug_rows = [
        r for r in rows
        if r["failure_class"] in {
            "assembly_shape_mismatch",
            "input_degree_incompatible",
            "enumeration_no_structures",
            "enumerator_degree_leak",
            "assembler_anchor_mismatch",
            "assembly_no_cif",
            "assembly_exception",
            "comparison_severe_difference",
            "comparison_warning_difference",
        }
    ]
    for row in debug_rows[:10]:
        labels = row.get("errors") or row.get("warnings") or row.get("diagnostic") or row.get("error_message")
        lines.append(f"- `{row['case']}`: `{row['failure_class']}`; {labels}")
    if not debug_rows:
        lines.append("- No generated-case debug target was found.")

    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark local reconstruction against ToBaCCo CIFs.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--limit", type=int, default=None, help="Run only the first N reference CIFs.")
    parser.add_argument("--case", action="append", default=[], help="Run a specific CIF file or case stem. Can repeat.")
    parser.add_argument("--max-cifs", type=int, default=1, help="Generated CIFs compared per reference case.")
    parser.add_argument("--max-results", type=int, default=5000, help="Enumeration result cap passed to main.run_enumeration.")
    args = parser.parse_args()

    out_dir = args.out_dir
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)

    ref_files = _scan_reference_files(limit=args.limit, cases=args.case)
    if not ref_files:
        raise SystemExit("No reference CIF files matched the benchmark selection.")

    start = time.time()
    rows = []
    for idx, ref_file in enumerate(ref_files, 1):
        print(f"[{idx}/{len(ref_files)}] {ref_file}")
        case, parse_error = _parse_case(ref_file)
        if parse_error:
            rows.append({
                "case": Path(ref_file).stem,
                "ref_file": ref_file,
                "topology": "",
                "nodes": "",
                "edge": "",
                "failure_class": "parse_error",
                "error_message": parse_error,
            })
            continue
        rows.append(_run_one_case(case, out_dir, args.max_cifs, args.max_results))

    elapsed = time.time() - start
    summary_rows = _topology_summary(rows)

    detail_csv = out_dir / "case_results.csv"
    summary_csv = out_dir / "topology_summary.csv"
    summary_md = out_dir / "benchmark_summary.md"
    detail_json = out_dir / "case_results.json"

    _write_csv(detail_csv, rows)
    _write_csv(summary_csv, summary_rows)
    detail_json.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n")
    _write_markdown(summary_md, rows, summary_rows, elapsed, args)

    print("\nBenchmark complete.")
    print(f"  Detail CSV: {detail_csv}")
    print(f"  Topology CSV: {summary_csv}")
    print(f"  Summary MD: {summary_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
