#!/usr/bin/env python3
"""Export a first ffc CIF batch for CO2/N2 adsorption screening."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

sys.path.insert(0, ".")

from assembler import StructureAssembler
from dft_validation.select_structures import classify_sym_classes, comp_label_key
from enumerator import MOFEnumerator
from main import calculate_topology_requirements, enrich_graph_data_with_topology, load_graph
from symmetry import SymmetryEngine


TOPOLOGY = "ffc"
FORCED_EDGE = "1B_4H_Ch"
DEFAULT_BATCH_ID = "batch_001"
DEFAULT_TARGET_COUNT = 200
DEFAULT_SEED = 42
MAX_ATOMS = 450

FIRST_BATCH_NODE_CHOICES = {
    3: ["3c_triazine_Ch", "3c_pyrrole_Ch", "3c_Cu_imadazoline_Ch", "3c_B_Ch"],
    4: [
        "4c_Cu_1_Ch",
        "4c_Cd_1_Ch",
        "4c_Ni_1_Ch",
        "4c_V_1_Ch",
        "4c_In_1_Ch",
        "4c_porphyrin_Ch",
        "4c_tetrathiafulvalene_Ch",
        "4c_B_Ch",
        "4c_pyrene_Ch",
        "4c_bicyclooctane_Ch",
    ],
}

EXCLUSION_RULES = {
    "3c_uranyl": "U-containing; force-field risk",
    "3c_uranyl_Ch": "U-containing; force-field risk",
    "4c_1anC_Ch": "Fr-containing placeholder-like chemistry",
    "4c_ethene_Ch": "actual anchor count is 2; not usable for ffc 4c nodes",
    "4c_Zn_tet": "X-only special anchor semantics; keep for manual review",
    "4c_pyrene_NU_Ch": "duplicate aromatic variant; omit from first batch",
}

DESCRIPTOR_ELEMENTS = ["H", "B", "C", "N", "O", "S", "Cu", "Cd", "Ni", "V", "In"]


@dataclass(frozen=True)
class NodeInfo:
    name: str
    degree_from_name: int | None
    inferred_connectivity: int
    path: str
    snapshot_name: str
    selected_for_batch: bool
    exclude_reason: str
    n_atoms: int
    n_heavy_atoms: int
    n_hetero_atoms: int
    n_metals: int
    element_counts: dict[str, int]
    family_tags: list[str]


def _project_path(*parts: str) -> Path:
    return Path(*parts)


def _degree_from_name(name: str) -> int | None:
    head = name.split("_", 1)[0]
    if head.endswith("c") and head[:-1].isdigit():
        return int(head[:-1])
    return None


def _read_atom_site(path: Path) -> list[tuple[str, str]]:
    lines = path.read_text(errors="replace").splitlines()
    atoms: list[tuple[str, str]] = []
    headers: list[str] = []
    in_atom_loop = False

    for line in lines:
        stripped = line.strip()
        if stripped == "loop_":
            if atoms:
                break
            headers = []
            in_atom_loop = False
            continue
        if stripped.startswith("_atom_site_"):
            headers.append(stripped.split()[0])
            in_atom_loop = True
            continue
        if in_atom_loop and stripped.startswith("_"):
            continue
        if not in_atom_loop or not stripped or stripped.startswith("#"):
            continue

        if "_atom_site_label" not in headers or "_atom_site_type_symbol" not in headers:
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        atoms.append((parts[0], parts[1]))

    return atoms


def _family_tags(name: str, element_counts: dict[str, int], n_atoms: int) -> list[str]:
    tags: list[str] = []
    lower = name.lower()
    if element_counts.get("N", 0) >= 3:
        tags.append("N-rich")
    elif element_counts.get("N", 0) > 0:
        tags.append("N-containing")
    if element_counts.get("O", 0) >= 4:
        tags.append("O-rich")
    if element_counts.get("S", 0) > 0:
        tags.append("S-containing")
    if any(el in element_counts for el in ("Cu", "Cd", "Ni", "V", "In", "Zn", "U", "Fr")):
        tags.append("metal-containing")
    if any(word in lower for word in ("triazine", "pyrrole", "imidazoline", "porphyrin")):
        tags.append("polar-ring")
    if any(word in lower for word in ("pyrene", "porphyrin", "fused")):
        tags.append("extended-aromatic")
    if n_atoms >= 20:
        tags.append("bulky")
    if not tags:
        tags.append("hydrocarbon")
    return sorted(set(tags))


def inspect_node_file(path: Path) -> NodeInfo:
    name = path.stem
    atoms = _read_atom_site(path)
    labels = [label for label, _ in atoms]
    element_counts = dict(Counter(symbol for _, symbol in atoms))
    inferred = sum(1 for label in labels if label.startswith("X"))
    degree = _degree_from_name(name)
    n_atoms = len(atoms)
    n_heavy = sum(count for element, count in element_counts.items() if element != "H")
    hetero = sum(element_counts.get(el, 0) for el in ("N", "O", "S", "F", "Cl", "Br", "I"))
    metals = sum(
        count
        for element, count in element_counts.items()
        if element not in {"H", "B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "X"}
    )

    selected = name in set(FIRST_BATCH_NODE_CHOICES.get(degree or -1, []))
    exclude_reason = EXCLUSION_RULES.get(name, "")
    if degree is None:
        selected = False
        exclude_reason = exclude_reason or "no connectivity prefix"
    elif inferred != degree:
        selected = False
        exclude_reason = exclude_reason or f"inferred connectivity {inferred} != filename degree {degree}"
    elif not selected and not exclude_reason:
        exclude_reason = "not in batch_001 node pool"

    snapshot_name = f"{name}_{inferred}c.cif"
    return NodeInfo(
        name=name,
        degree_from_name=degree,
        inferred_connectivity=inferred,
        path=str(path),
        snapshot_name=snapshot_name,
        selected_for_batch=selected,
        exclude_reason=exclude_reason,
        n_atoms=n_atoms,
        n_heavy_atoms=n_heavy,
        n_hetero_atoms=hetero,
        n_metals=metals,
        element_counts=element_counts,
        family_tags=_family_tags(name, element_counts, n_atoms),
    )


def screen_nodes(nodes_dir: Path = Path("nodes")) -> list[NodeInfo]:
    paths = sorted(nodes_dir.glob("3c*.cif")) + sorted(nodes_dir.glob("4c*.cif"))
    return [inspect_node_file(path) for path in paths]


def selected_atom_choices(node_infos: Iterable[NodeInfo]) -> dict[int, list[str]]:
    choices: dict[int, list[str]] = defaultdict(list)
    for info in node_infos:
        if info.selected_for_batch and info.degree_from_name in (3, 4):
            choices[int(info.degree_from_name)].append(info.name)
    return {degree: sorted(names) for degree, names in sorted(choices.items())}


def _node_info_map(node_infos: Iterable[NodeInfo]) -> dict[str, NodeInfo]:
    return {info.name: info for info in node_infos}


def _composition_dict(coloring: list[str]) -> dict[str, int]:
    return dict(sorted(Counter(coloring).items()))


def _composition_string(coloring: list[str]) -> str:
    return comp_label_key(coloring)


def _candidate_tags(coloring: list[str], node_map: dict[str, NodeInfo]) -> list[str]:
    tags: set[str] = set()
    for name in coloring:
        tags.update(node_map[name].family_tags)
    return sorted(tags)


def _candidate_element_counts(coloring: list[str], node_map: dict[str, NodeInfo]) -> Counter:
    counts: Counter = Counter()
    for name in coloring:
        counts.update(node_map[name].element_counts)
    return counts


def _descriptor_vector(record: dict, node_names: list[str], node_map: dict[str, NodeInfo]) -> list[float]:
    coloring = record["coloring"]
    composition = Counter(coloring)
    elements = _candidate_element_counts(coloring, node_map)
    tags = set(_candidate_tags(coloring, node_map))
    n_atoms = sum(node_map[name].n_atoms for name in coloring)
    hetero = sum(node_map[name].n_hetero_atoms for name in coloring)
    metals = sum(node_map[name].n_metals for name in coloring)

    vector: list[float] = [float(composition.get(name, 0)) for name in node_names]
    vector.extend(float(elements.get(el, 0)) for el in DESCRIPTOR_ELEMENTS)
    vector.extend(
        [
            float(n_atoms),
            float(hetero),
            float(metals),
            float(record["rho"]),
            float("N-rich" in tags or "N-containing" in tags),
            float("O-rich" in tags),
            float("S-containing" in tags),
            float("extended-aromatic" in tags),
            float("bulky" in tags),
        ]
    )
    return vector


def _normalize_vectors(vectors: list[list[float]]) -> list[list[float]]:
    if not vectors:
        return []
    dims = len(vectors[0])
    mins = [min(v[i] for v in vectors) for i in range(dims)]
    maxs = [max(v[i] for v in vectors) for i in range(dims)]
    out: list[list[float]] = []
    for vector in vectors:
        norm = []
        for i, value in enumerate(vector):
            span = maxs[i] - mins[i]
            norm.append(0.0 if span == 0 else (value - mins[i]) / span)
        out.append(norm)
    return out


def _distance(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _enumerate_ffc(atom_choices: dict[int, list[str]]) -> tuple[list[dict], dict, int]:
    graph_data = enrich_graph_data_with_topology(load_graph("data/ffc/ffc_graph.json"))
    topo_requirements = calculate_topology_requirements(graph_data)
    node_degrees = [node.get("degree", 0) for node in graph_data["nodes"]]
    perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
    ligands = sorted({ligand for names in atom_choices.values() for ligand in names})
    results = MOFEnumerator(
        perms,
        graph_data["n_nodes"],
        verbose=False,
        topology_requirements=topo_requirements,
        user_selected_ligands=atom_choices,
        node_degrees=node_degrees,
    ).enumerate(ligands)

    gamma = len(perms)
    classified = classify_sym_classes(results)
    records = []
    for result in classified:
        records.append(
            {
                "source_enum_id": result["id"],
                "topology": TOPOLOGY,
                "composition": _composition_string(result["coloring"]),
                "composition_dict": _composition_dict(result["coloring"]),
                "coloring": result["coloring"],
                "stabilizer_size": result["stabilizer_size"],
                "gamma_size": gamma,
                "rho": result["stabilizer_size"] / gamma,
                "sym_class": result["sym_class"],
            }
        )
    return records, graph_data, gamma


def _sorted_bucket(records: Iterable[dict]) -> list[dict]:
    return sorted(
        records,
        key=lambda rec: (
            rec["composition"],
            -rec["rho"],
            rec["source_enum_id"],
            ",".join(rec["coloring"]),
        ),
    )


def _balanced_take(records: list[dict], count: int) -> list[dict]:
    by_comp: dict[str, list[dict]] = defaultdict(list)
    for record in _sorted_bucket(records):
        by_comp[record["composition"]].append(record)

    selected: list[dict] = []
    used_ids: set[int] = set()
    comp_keys = sorted(by_comp)
    while len(selected) < count:
        made_progress = False
        for key in comp_keys:
            bucket = by_comp[key]
            while bucket and id(bucket[0]) in used_ids:
                bucket.pop(0)
            if not bucket:
                continue
            record = bucket.pop(0)
            selected.append(record)
            used_ids.add(id(record))
            made_progress = True
            if len(selected) == count:
                break
        if not made_progress:
            break
    return selected


def _max_min_diverse(records: list[dict], count: int, node_names: list[str],
                     node_map: dict[str, NodeInfo], seed_records: list[dict]) -> list[dict]:
    if not records or count <= 0:
        return []

    candidates = _sorted_bucket(records)
    all_records = seed_records + candidates
    vectors = [_descriptor_vector(record, node_names, node_map) for record in all_records]
    norm = np.array(_normalize_vectors(vectors), dtype=float)
    seed_vectors = norm[:len(seed_records)]
    candidate_vectors = norm[len(seed_records):]

    chosen: list[dict] = []
    available = np.ones(len(candidates), dtype=bool)
    min_dist = np.full(len(candidates), np.inf, dtype=float)

    for ref in seed_vectors:
        dist = np.linalg.norm(candidate_vectors - ref, axis=1)
        min_dist = np.minimum(min_dist, dist)

    if len(seed_vectors) == 0 and len(candidates) > 0:
        first = max(
            range(len(candidates)),
            key=lambda i: (
                candidates[i]["rho"],
                sum(candidate_vectors[i]),
                -candidates[i]["source_enum_id"],
            ),
        )
        chosen.append(candidates[first])
        available[first] = False
        dist = np.linalg.norm(candidate_vectors - candidate_vectors[first], axis=1)
        min_dist = np.minimum(min_dist, dist)

    while len(chosen) < count and available.any():
        valid = np.flatnonzero(available)
        max_dist = min_dist[valid].max()
        tied = valid[np.isclose(min_dist[valid], max_dist)]
        best_idx = max(
            tied.tolist(),
            key=lambda i: (candidates[i]["rho"], -candidates[i]["source_enum_id"]),
        )
        chosen.append(candidates[best_idx])
        available[best_idx] = False
        dist = np.linalg.norm(candidate_vectors - candidate_vectors[best_idx], axis=1)
        min_dist = np.minimum(min_dist, dist)

    return chosen


def build_selection_plan(target_count: int = DEFAULT_TARGET_COUNT,
                         seed: int = DEFAULT_SEED,
                         include_backups: bool = False) -> tuple[list[dict], list[NodeInfo], dict, int]:
    if seed != DEFAULT_SEED:
        raise ValueError("Only seed=42 is currently supported for reproducible batch_001 ordering")
    if target_count != DEFAULT_TARGET_COUNT:
        raise ValueError("batch_001 selection ratios are fixed for target_count=200")

    node_infos = screen_nodes()
    atom_choices = selected_atom_choices(node_infos)
    records, graph_data, gamma = _enumerate_ffc(atom_choices)
    node_map = _node_info_map(node_infos)
    node_names = sorted({name for names in atom_choices.values() for name in names})

    high = [record for record in records if record["sym_class"] == "High"]
    med = [record for record in records if record["sym_class"] == "Med"]
    low = [record for record in records if record["sym_class"] == "Low"]

    high_quota = 80 if not include_backups else 110
    med_quota = 60 if not include_backups else 90
    diversity_quota = 40 if not include_backups else 90
    low_quota = 20 if not include_backups else 60

    symmetry_high = _balanced_take(high, high_quota)
    selected_keys = {record_key(record) for record in symmetry_high}
    symmetry_med = _balanced_take(
        [record for record in med if record_key(record) not in selected_keys],
        med_quota,
    )
    symmetry_order = (
        symmetry_high[:80]
        + symmetry_med[:60]
        + symmetry_high[80:]
        + symmetry_med[60:]
    )
    symmetry_records = [
        {**record, "selection_bucket": "symmetry_priority"}
        for record in symmetry_order
    ]

    selected_keys = {record_key(record) for record in symmetry_records}
    diversity_pool = [
        record for record in records
        if record["sym_class"] != "Low" and record_key(record) not in selected_keys
    ]
    diversity_raw = _max_min_diverse(diversity_pool, diversity_quota, node_names, node_map, symmetry_records)
    if len(diversity_raw) < diversity_quota:
        fallback_pool = [
            record for record in records
            if record_key(record) not in selected_keys
            and record_key(record) not in {record_key(r) for r in diversity_raw}
        ]
        diversity_raw.extend(
            _max_min_diverse(fallback_pool, diversity_quota - len(diversity_raw), node_names, node_map,
                             symmetry_records + diversity_raw)
        )
    diversity_records = [
        {**record, "selection_bucket": "descriptor_diversity"}
        for record in diversity_raw[:diversity_quota]
    ]

    selected_keys = {record_key(record) for record in symmetry_records + diversity_records}
    low_records = [
        {**record, "selection_bucket": "low_sym_control"}
        for record in _balanced_take(
            [record for record in low if record_key(record) not in selected_keys],
            low_quota,
        )
    ]

    plan = symmetry_records + diversity_records + low_records
    expected = target_count if not include_backups else high_quota + med_quota + diversity_quota + low_quota
    if len(plan) != expected:
        raise RuntimeError(f"Selection produced {len(plan)} records, expected {expected}")
    return plan, node_infos, graph_data, gamma


def record_key(record: dict) -> tuple:
    return (
        record["source_enum_id"],
        tuple(record["coloring"]),
        record["stabilizer_size"],
    )


def _copy_node_snapshots(node_infos: list[NodeInfo], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for info in node_infos:
        target = output_dir / info.snapshot_name
        shutil.copyfile(info.path, target)


def _node_rows(node_infos: list[NodeInfo]) -> list[dict]:
    rows = []
    for info in node_infos:
        rows.append(
            {
                "node_name": info.name,
                "degree_from_name": info.degree_from_name,
                "inferred_connectivity": info.inferred_connectivity,
                "snapshot_name": info.snapshot_name,
                "selected_for_batch": info.selected_for_batch,
                "exclude_reason": info.exclude_reason,
                "n_atoms": info.n_atoms,
                "n_heavy_atoms": info.n_heavy_atoms,
                "n_hetero_atoms": info.n_hetero_atoms,
                "n_metals": info.n_metals,
                "element_counts_json": json.dumps(info.element_counts, sort_keys=True),
                "family_tags": ";".join(info.family_tags),
                "source_path": info.path,
            }
        )
    return rows


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


VALIDATION_FIELD_DEFAULTS = {
    "validation_passed": "",
    "clash_passed": "",
    "topology_passed": "",
    "symmetry_passed": "",
    "n_clashes": "",
    "worst_clash_ratio": "",
    "clash_message": "",
    "topology_message": "",
    "symmetry_message": "",
    "detected_sg": "",
    "validation_warning": "",
    "validation_json": "",
}


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _check_result(validation: dict, check_name: str) -> dict:
    result = validation.get(check_name, {})
    return result if isinstance(result, dict) else {}


def _check_passed(validation: dict, check_name: str) -> bool | str:
    result = _check_result(validation, check_name)
    if result.get("skipped"):
        return ""
    passed = result.get("passed")
    return "" if passed is None else bool(passed)


def _failed_validation_checks(validation: dict) -> list[str]:
    failed = []
    for check_name in ("clash", "topology", "symmetry"):
        result = _check_result(validation, check_name)
        if result.get("skipped"):
            continue
        if result.get("passed") is False:
            failed.append(check_name)
    return failed


def _validation_fields(metadata: dict | None) -> dict:
    if not metadata:
        return dict(VALIDATION_FIELD_DEFAULTS)

    validation = metadata.get("validation", {})
    validation = validation if isinstance(validation, dict) else {}
    clash = _check_result(validation, "clash")
    topology = _check_result(validation, "topology")
    symmetry = _check_result(validation, "symmetry")
    worst_clash = clash.get("worst_clash", {})
    if not isinstance(worst_clash, dict):
        worst_clash = {}
    failed = _failed_validation_checks(validation)

    return {
        "validation_passed": bool(metadata.get("validation_passed", validation.get("passed", False))),
        "clash_passed": _check_passed(validation, "clash"),
        "topology_passed": _check_passed(validation, "topology"),
        "symmetry_passed": _check_passed(validation, "symmetry"),
        "n_clashes": clash.get("n_clashes", ""),
        "worst_clash_ratio": worst_clash.get("ratio", ""),
        "clash_message": clash.get("message", ""),
        "topology_message": topology.get("message", ""),
        "symmetry_message": symmetry.get("message", ""),
        "detected_sg": metadata.get("detected_sg", symmetry.get("detected_sg", "")),
        "validation_warning": ";".join(failed),
        "validation_json": json.dumps(validation, ensure_ascii=False, sort_keys=True, default=_json_default),
    }


def _manifest_row(record: dict, candidate_id: str, cif_path: Path, structure,
                  node_map: dict[str, NodeInfo], metadata: dict | None = None) -> dict:
    counts = Counter(record["coloring"])
    source_node_files = [node_map[name].path for name in sorted(counts)]
    row = {
        "candidate_id": candidate_id,
        "topology": TOPOLOGY,
        "composition": record["composition"],
        "rho": f"{record['rho']:.6f}",
        "sym_class": record["sym_class"],
        "cif_path": str(cif_path),
        "selection_bucket": record["selection_bucket"],
        "source_enum_id": record["source_enum_id"],
        "stabilizer_size": record["stabilizer_size"],
        "gamma_size": record["gamma_size"],
        "coloring_json": json.dumps(record["coloring"], sort_keys=False),
        "n_atoms": len(structure),
        "formula": structure.composition.reduced_formula,
        "density_g_cm3": f"{structure.density:.6f}",
        "node_family_tags": ";".join(_candidate_tags(record["coloring"], node_map)),
        "source_node_files": json.dumps(source_node_files, sort_keys=True),
    }
    row.update(_validation_fields(metadata))
    return row


def _assemble_manifest(plan: list[dict], graph_data: dict, node_infos: list[NodeInfo],
                       out_dir: Path, target_count: int) -> tuple[list[dict], list[dict]]:
    from pymatgen.core import Structure

    cif_dir = out_dir / "cifs"
    cif_dir.mkdir(parents=True, exist_ok=True)
    assembler = StructureAssembler(output_dir=str(cif_dir))
    node_map = _node_info_map(node_infos)
    quotas = {
        "symmetry_priority": 140,
        "descriptor_diversity": 40,
        "low_sym_control": 20,
    }
    successes_by_bucket = Counter()
    manifest: list[dict] = []
    failures: list[dict] = []

    for record in plan:
        bucket = record["selection_bucket"]
        if successes_by_bucket[bucket] >= quotas[bucket]:
            continue
        candidate_number = len(manifest) + 1
        candidate_id = f"ffc_co2n2_batch001_{candidate_number:04d}"
        cif_path = cif_dir / f"{candidate_id}_{bucket}_rho{record['rho']:.4f}_{record['sym_class']}.cif"
        structure_map = {str(i): name for i, name in enumerate(record["coloring"])}
        try:
            metadata = assembler.assemble_and_save(
                str(cif_path),
                structure_map,
                graph_data,
                forced_edge=FORCED_EDGE,
                return_metadata=True,
            )
            structure = Structure.from_file(cif_path)
            if len(structure) > MAX_ATOMS:
                failures.append(
                    {
                        "source_enum_id": record["source_enum_id"],
                        "selection_bucket": bucket,
                        "reason": f"n_atoms {len(structure)} > {MAX_ATOMS}",
                    }
                )
                cif_path.unlink(missing_ok=True)
                continue
        except Exception as exc:
            failures.append(
                {
                    "source_enum_id": record["source_enum_id"],
                    "selection_bucket": bucket,
                    "reason": str(exc),
                }
            )
            continue

        manifest.append(_manifest_row(record, candidate_id, cif_path, structure, node_map, metadata))
        successes_by_bucket[bucket] += 1
        if len(manifest) >= target_count:
            break

    if len(manifest) != target_count:
        raise RuntimeError(
            f"Assembled {len(manifest)} candidates, expected {target_count}. "
            f"Failures: {len(failures)}"
        )
    return manifest, failures


def _write_summary(out_dir: Path, manifest: list[dict], failures: list[dict],
                   node_infos: list[NodeInfo], gamma: int) -> None:
    by_bucket = Counter(row["selection_bucket"] for row in manifest)
    by_class = Counter(row["sym_class"] for row in manifest)
    selected_nodes = [info.name for info in node_infos if info.selected_for_batch]
    summary = {
        "batch_id": out_dir.name,
        "topology": TOPOLOGY,
        "target_count": len(manifest),
        "gamma_size": gamma,
        "max_atoms": MAX_ATOMS,
        "bucket_counts": dict(sorted(by_bucket.items())),
        "sym_class_counts": dict(sorted(by_class.items())),
        "selected_nodes": selected_nodes,
        "excluded_nodes": {
            info.name: info.exclude_reason
            for info in node_infos
            if info.exclude_reason and not info.selected_for_batch
        },
        "failed_assembly_count": len(failures),
        "failed_assemblies": failures,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    lines = [
        f"# {out_dir.name} ffc CO2/N2 candidate export",
        "",
        f"- Topology: `{TOPOLOGY}`",
        f"- Exported CIFs: `{len(manifest)}`",
        f"- Automorphism group size: `{gamma}`",
        f"- Max atom cap: `{MAX_ATOMS}`",
        f"- Bucket counts: `{dict(sorted(by_bucket.items()))}`",
        f"- Symmetry classes: `{dict(sorted(by_class.items()))}`",
        f"- Failed / skipped assemblies: `{len(failures)}`",
        "",
        "This batch prepares CIFs and metadata for later CO2/N2 screening; it does not run GCMC.",
    ]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")


def export_batch(batch_id: str = DEFAULT_BATCH_ID, target_count: int = DEFAULT_TARGET_COUNT,
                 out_root: Path = Path("results/co2n2_ffc_screening"),
                 dry_run: bool = False) -> list[dict]:
    plan, node_infos, graph_data, gamma = build_selection_plan(
        target_count=target_count,
        include_backups=not dry_run,
    )
    if dry_run:
        rows = []
        for idx, record in enumerate(plan, start=1):
            rows.append(
                {
                    "candidate_id": f"ffc_co2n2_batch001_{idx:04d}",
                    "topology": TOPOLOGY,
                    "composition": record["composition"],
                    "rho": f"{record['rho']:.6f}",
                    "sym_class": record["sym_class"],
                    "cif_path": f"{out_root}/{batch_id}/cifs/ffc_co2n2_batch001_{idx:04d}.cif",
                    "selection_bucket": record["selection_bucket"],
                    "source_enum_id": record["source_enum_id"],
                    "stabilizer_size": record["stabilizer_size"],
                    "gamma_size": record["gamma_size"],
                    "coloring_json": json.dumps(record["coloring"]),
                    "n_atoms": "",
                    "formula": "",
                    "density_g_cm3": "",
                    "node_family_tags": "",
                    "source_node_files": "",
                    **VALIDATION_FIELD_DEFAULTS,
                }
            )
        return rows

    out_dir = out_root / batch_id
    out_dir.mkdir(parents=True, exist_ok=True)
    _copy_node_snapshots(node_infos, out_dir / "nodes_snapshot")
    _write_csv(out_dir / "nodes_screening.csv", _node_rows(node_infos))
    manifest, failures = _assemble_manifest(plan, graph_data, node_infos, out_dir, target_count)
    _write_csv(out_dir / "manifest.csv", manifest)
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    _write_summary(out_dir, manifest, failures, node_infos, gamma)
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-id", default=DEFAULT_BATCH_ID)
    parser.add_argument("--target-count", type=int, default=DEFAULT_TARGET_COUNT)
    parser.add_argument("--out-root", default="results/co2n2_ffc_screening")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    rows = export_batch(
        batch_id=args.batch_id,
        target_count=args.target_count,
        out_root=Path(args.out_root),
        dry_run=args.dry_run,
    )
    print(f"[co2n2-export] prepared {len(rows)} candidates")
    if not args.dry_run:
        print(f"[co2n2-export] output: {Path(args.out_root) / args.batch_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
