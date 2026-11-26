"""Audit whether ToBaCCo reference CIFs match local topology graph semantics.

The benchmark filename tells us which building blocks ToBaCCo used.  This
script checks that claim against the reference CIF composition and then compares
the inferred node/edge counts with the local ``data/<topology>/*_graph.json``.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from itertools import product
from pathlib import Path

from pymatgen.core import Structure

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from main import calculate_topology_requirements, load_graph  # noqa: E402
from tests.name_parser import parse_filename  # noqa: E402


REF_DIR = PROJECT_ROOT / "tests" / "tobacco_output_cifs"
NODE_DIR = PROJECT_ROOT / "nodes"
EDGE_DIR = PROJECT_ROOT / "edges"
DATA_DIR = PROJECT_ROOT / "data"


def _json_dump(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _declared_degree(block_name: str) -> int | None:
    match = re.match(r"(\d+)c_", block_name)
    return int(match.group(1)) if match else None


def _read_cif_atom_counts(path: Path, exclude_x: bool = True) -> dict[str, int]:
    """Read atom counts from the atom_site loop, optionally ignoring X anchors."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    counts: Counter[str] = Counter()
    in_loop = False
    headers: list[str] = []
    label_idx = type_idx = None

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line == "loop_":
            in_loop = True
            headers = []
            label_idx = type_idx = None
            continue
        if not in_loop:
            continue
        if line.startswith("_"):
            headers.append(line)
            if line == "_atom_site_label":
                label_idx = len(headers) - 1
            elif line == "_atom_site_type_symbol":
                type_idx = len(headers) - 1
            continue
        if label_idx is None or type_idx is None:
            if headers:
                in_loop = False
            continue
        parts = line.split()
        if len(parts) <= max(label_idx, type_idx):
            continue
        label = parts[label_idx]
        element = parts[type_idx]
        if exclude_x and label.startswith("X"):
            continue
        counts[element] += 1

    return dict(counts)


def _pymatgen_counts(path: Path) -> dict[str, int]:
    structure = Structure.from_file(str(path))
    return {el: int(round(amount)) for el, amount in structure.composition.as_dict().items()}


def _block_composition(block_name: str, block_type: str) -> dict[str, int]:
    base_dir = NODE_DIR if block_type == "node" else EDGE_DIR
    path = base_dir / f"{block_name}.cif"
    if not path.exists():
        raise FileNotFoundError(path)
    # ToBaCCo's reference CIFs retain the atoms labelled X as real C sites in
    # the final composition, so composition-based block inference must include
    # them.  Connectivity is still interpreted separately from the filename and
    # the number of X anchors in the block CIF.
    return _read_cif_atom_counts(path, exclude_x=False)


def _composition_bounds(ref_counts: dict[str, int], variables: list[str],
                        formulas: dict[str, dict[str, int]]) -> list[range]:
    bounds = []
    total_atoms = sum(ref_counts.values())
    for var in variables:
        formula = formulas[var]
        candidates = [
            ref_counts.get(el, 0) // count
            for el, count in formula.items()
            if count > 0
        ]
        upper = min(candidates) if candidates else total_atoms
        bounds.append(range(upper + 1))
    return bounds


def _infer_block_counts(ref_counts: dict[str, int], node_names: list[str],
                        edge_name: str) -> tuple[dict[str, int], int | None, str]:
    variables = sorted(set(node_names)) + [edge_name]
    formulas: dict[str, dict[str, int]] = {}
    for node_name in sorted(set(node_names)):
        formulas[node_name] = _block_composition(node_name, "node")
    formulas[edge_name] = _block_composition(edge_name, "edge")

    solutions = []
    for combo in product(*_composition_bounds(ref_counts, variables, formulas)):
        total: Counter[str] = Counter()
        for var, count in zip(variables, combo):
            for el, amount in formulas[var].items():
                total[el] += amount * count
        if dict(total) == ref_counts:
            solutions.append(dict(zip(variables, combo)))
            if len(solutions) > 1:
                break

    if not solutions:
        return {}, None, "unsolved"
    if len(solutions) > 1:
        return solutions[0], solutions[0].get(edge_name), "ambiguous"
    solution = solutions[0]
    return solution, solution.get(edge_name), "unique"


def _local_graph_summary(topology: str) -> tuple[dict[str, object], str]:
    graph_path = DATA_DIR / topology / f"{topology}_graph.json"
    if not graph_path.exists():
        return {}, "missing_topology_graph"

    graph_data = load_graph(str(graph_path))
    requirements = calculate_topology_requirements(graph_data)
    return {
        "local_nodes": len(graph_data["nodes"]),
        "local_edges": len(graph_data.get("edges", [])),
        "local_degree_requirements": requirements,
    }, ""


def _semantic_class(local: dict[str, object], inferred: dict[str, int],
                    edge_name: str, inference_status: str) -> str:
    if not local:
        return "missing_topology_graph"
    if inference_status == "unsolved":
        return "composition_inference_failed"

    inferred_node_counts = {
        _declared_degree(name): count
        for name, count in inferred.items()
        if name != edge_name
    }
    inferred_node_counts = {
        int(deg): count for deg, count in inferred_node_counts.items()
        if deg is not None
    }
    local_req = {int(k): int(v) for k, v in local["local_degree_requirements"].items()}
    inferred_nodes_total = sum(
        count for name, count in inferred.items() if name != edge_name
    )
    inferred_edges = inferred.get(edge_name)

    checks = [
        int(local["local_nodes"]) == inferred_nodes_total,
        int(local["local_edges"]) == inferred_edges,
        local_req == inferred_node_counts,
    ]
    return "semantic_match" if all(checks) else "semantic_mismatch"


def audit_case(ref_path: Path) -> dict[str, str]:
    topology, node_names, edge_name = parse_filename(ref_path.name)
    local, input_issue = _local_graph_summary(topology)
    ref_counts = _pymatgen_counts(ref_path)
    inferred, inferred_edges, inference_status = _infer_block_counts(
        ref_counts, node_names, edge_name
    )

    inferred_node_counts = {
        name: count for name, count in inferred.items() if name != edge_name
    }
    inferred_degree_requirements: Counter[int] = Counter()
    for name, count in inferred_node_counts.items():
        degree = _declared_degree(name)
        if degree is not None:
            inferred_degree_requirements[degree] += count

    semantic_class = _semantic_class(local, inferred, edge_name, inference_status)
    row = {
        "ref_file": ref_path.name,
        "topology": topology,
        "filename_nodes": ";".join(node_names),
        "filename_edge": edge_name,
        "input_issue": input_issue,
        "ref_composition": _json_dump(ref_counts),
        "inference_status": inference_status,
        "inferred_node_counts": _json_dump(inferred_node_counts),
        "inferred_edge_count": "" if inferred_edges is None else str(inferred_edges),
        "inferred_degree_requirements": _json_dump(dict(sorted(inferred_degree_requirements.items()))),
        "local_nodes": str(local.get("local_nodes", "")),
        "local_edges": str(local.get("local_edges", "")),
        "local_degree_requirements": _json_dump(local.get("local_degree_requirements", {})),
        "semantic_class": semantic_class,
    }
    return row


def write_outputs(rows: list[dict[str, str]], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else []

    csv_path = out_dir / "topology_semantics_audit.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    json_path = out_dir / "topology_semantics_audit.json"
    json_path.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")

    counts = Counter(row["semantic_class"] for row in rows)
    md_lines = [
        "# ToBaCCo / Local Topology Semantics Audit",
        "",
        f"- Cases: `{len(rows)}`",
        f"- Semantic matches: `{counts.get('semantic_match', 0)}`",
        f"- Semantic mismatches: `{counts.get('semantic_mismatch', 0)}`",
        f"- Missing local topology graphs: `{counts.get('missing_topology_graph', 0)}`",
        f"- Composition inference failures: `{counts.get('composition_inference_failed', 0)}`",
        "",
        "| Topology | Reference | Local n/e | Inferred nodes | Inferred edges | Local degree req | Inferred degree req | Class |",
        "| --- | --- | ---: | --- | ---: | --- | --- | --- |",
    ]
    for row in rows:
        local_ne = f"{row['local_nodes']}/{row['local_edges']}" if row["local_nodes"] else ""
        md_lines.append(
            f"| {row['topology']} | `{row['ref_file']}` | {local_ne} | "
            f"`{row['inferred_node_counts']}` | {row['inferred_edge_count']} | "
            f"`{row['local_degree_requirements']}` | "
            f"`{row['inferred_degree_requirements']}` | {row['semantic_class']} |"
        )

    (out_dir / "topology_semantics_audit.md").write_text(
        "\n".join(md_lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-dir", type=Path, default=REF_DIR)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "tobacco_topology_semantics_audit",
    )
    args = parser.parse_args()

    rows = [audit_case(path) for path in sorted(args.ref_dir.glob("*.cif"))]
    write_outputs(rows, args.out_dir)

    counts = Counter(row["semantic_class"] for row in rows)
    print(f"Wrote {len(rows)} rows to {args.out_dir}")
    print(dict(sorted(counts.items())))


if __name__ == "__main__":
    main()
