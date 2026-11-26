import argparse
import csv
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
from pymatgen.core import Structure

sys.path.insert(0, ".")

from dft_validation import config
from dft_validation.cp2k_runner import (
    generate_cp2k_input,
    parse_energy_hartrees,
    parse_relaxation_completion,
    run_one,
    write_relaxed_cif_from_restart,
)
from structure_validator import _get_covalent_radius


FIG2_CASES = [
    {
        "id": "fig2_ffc_clash_flagged",
        "topology": "ffc",
        "source_role": "fig2_or_reviewer_clash",
        "rho": 0.5,
        "sym_class": "High",
        "source_cif": "results/reconstruction_summary_split/ffc/ffc/ffc_High_rho0.5000_004.cif",
    },
    {
        "id": "fig2_fjh_representative",
        "topology": "fjh",
        "source_role": "fig2_representative",
        "rho": 1.0,
        "sym_class": "High",
        "source_cif": "results/reconstruction_summary_split/fjh/fjh/fjh_High_rho1.0000_000.cif",
    },
    {
        "id": "fig2_vba_representative",
        "topology": "vba",
        "source_role": "fig2_representative",
        "rho": 0.5,
        "sym_class": "High",
        "source_cif": "results/reconstruction_summary_split/vba/vba/vba_High_rho0.5000_000.cif",
    },
    {
        "id": "fig2_hwx_representative",
        "topology": "hwx",
        "source_role": "fig2_representative",
        "rho": 1.0,
        "sym_class": "High",
        "source_cif": "results/reconstruction_summary_split/hwx/hwx/hwx_High_rho1.0000_000.cif",
    },
    {
        "id": "fig2_nbo_stress_case",
        "topology": "nbo",
        "source_role": "fig2_or_reviewer_stress",
        "rho": 1.0,
        "sym_class": "High",
        "source_cif": "results/reconstruction_summary_split/nbo/nbo/nbo_High_rho1.0000_000.cif",
    },
]

DFT_MINIMUM_IDS = [
    "ffc_comp_a_11",
    "ffc_comp_b_21",
    "fjh_comp_a_12",
    "fjh_comp_b_1",
]


def _load_single_point_records() -> dict:
    path = Path(config.RESULTS_DIR) / "selected_structures.json"
    if not path.exists():
        return {}
    return {rec["id"]: rec for rec in json.loads(path.read_text())}


def _copy_initial_cif(src: str, rec_id: str, results_dir: Path) -> str:
    initial_dir = results_dir / "initial"
    initial_dir.mkdir(parents=True, exist_ok=True)
    dst = initial_dir / f"{rec_id}.cif"
    shutil.copy(src, dst)
    return str(dst)


def _record_from_cif(case: dict, results_dir: Path) -> dict:
    src = Path(case["source_cif"])
    if not src.exists():
        raise FileNotFoundError(f"missing CIF for {case['id']}: {src}")
    struct = Structure.from_file(str(src))
    initial_cif = _copy_initial_cif(str(src), case["id"], results_dir)
    elems = sorted({str(site.specie.symbol) for site in struct})
    return {
        "id": case["id"],
        "topology": case["topology"],
        "source_role": case["source_role"],
        "rho": case.get("rho"),
        "sym_class": case.get("sym_class"),
        "init_cif": initial_cif,
        "clean_cif": initial_cif,
        "source_cif": str(src),
        "n_atoms": len(struct),
        "formula": struct.composition.reduced_formula,
        "elements": elems,
        "status": "selected",
    }


def _record_from_single_point(rec_id: str, source_records: dict, results_dir: Path) -> dict:
    old = source_records.get(rec_id)
    if old is None:
        raise KeyError(f"{rec_id} not found in {config.RESULTS_DIR}/selected_structures.json")
    src = old.get("clean_cif") or old.get("init_cif")
    if not src or not Path(src).exists():
        raise FileNotFoundError(f"missing clean/init CIF for {rec_id}: {src}")
    initial_cif = _copy_initial_cif(src, rec_id, results_dir)
    struct = Structure.from_file(initial_cif)
    return {
        "id": rec_id,
        "topology": old.get("topology"),
        "composition_label": old.get("composition_label"),
        "source_role": "single_point_minimum",
        "rho": old.get("rho"),
        "sym_class": old.get("sym_class"),
        "init_cif": initial_cif,
        "clean_cif": initial_cif,
        "source_cif": src,
        "n_atoms": len(struct),
        "formula": struct.composition.reduced_formula,
        "elements": sorted({str(site.specie.symbol) for site in struct}),
        "single_point_e_per_atom": old.get("dft_e_per_atom"),
        "single_point_e_hartree": old.get("dft_e_hartree"),
        "status": "selected",
    }


def prepare_records(preset: str, results_dir: Path, only_id: str | None = None) -> list[dict]:
    source_records = _load_single_point_records()
    records = []
    include_all = only_id is None

    if preset in {"fig2", "reviewer7", "reviewer9"}:
        fig2_cases = [
            case for case in FIG2_CASES
            if include_all or case["id"] == only_id
        ]
        records.extend(_record_from_cif(case, results_dir) for case in fig2_cases)

    if preset == "ffc-smoke":
        ids = ["ffc_comp_a_11", "ffc_comp_b_21"]
    elif preset == "reviewer7":
        ids = ["ffc_comp_a_11", "ffc_comp_b_21"]
    elif preset == "reviewer9":
        ids = DFT_MINIMUM_IDS
    elif preset == "fig2":
        ids = []
    else:
        raise ValueError(f"unknown preset: {preset}")

    if only_id is not None:
        ids = [rec_id for rec_id in ids if rec_id == only_id]
        if not records and not ids:
            raise ValueError(f"--only {only_id!r} is not part of preset {preset!r}")

    records.extend(_record_from_single_point(x, source_records, results_dir) for x in ids)
    return records


def _selection_path(results_dir: Path) -> Path:
    return results_dir / "relaxation_selection.json"


def save_records(records: list[dict], results_dir: Path) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    _selection_path(results_dir).write_text(json.dumps(records, indent=2))
    write_manifest(records, results_dir)


def load_records(results_dir: Path) -> list[dict]:
    path = _selection_path(results_dir)
    if not path.exists():
        raise FileNotFoundError(f"missing relaxation selection: {path}")
    return json.loads(path.read_text())


def write_manifest(records: list[dict], results_dir: Path) -> None:
    fields = [
        "id", "topology", "source_role", "composition_label", "rho", "sym_class",
        "n_atoms", "formula", "elements", "status", "source_cif", "clean_cif",
        "relaxed_cif", "relax_run_type", "relaxed_extract_status",
        "cp2k_program_ended", "relaxation_completed_marker",
        "optimization_steps", "relaxation_convergence_flags_all_yes",
        "relaxation_convergence_flags",
        "initial_e_hartree", "final_e_hartree", "delta_e_hartree",
        "delta_e_ev_per_atom", "volume_change_pct", "rmsd_a",
        "bond_signature_retained",
    ]
    path = results_dir / "relaxation_manifest.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for rec in records:
            row = {k: rec.get(k) for k in fields}
            row["elements"] = ";".join(rec.get("elements") or [])
            writer.writerow(row)


def write_inputs(records: list[dict], results_dir: Path, run_type: str, cutoff_ry: int | None) -> list[dict]:
    updated = []
    for rec in records:
        run_dir = results_dir / rec["topology"] / "runs" / rec["id"]
        run_dir.mkdir(parents=True, exist_ok=True)
        try:
            inp = generate_cp2k_input(
                rec, str(run_dir), rec["topology"],
                relax_mode_override=run_type, cutoff_ry=cutoff_ry,
            )
            updated.append({
                **rec,
                "relax_run_type": run_type,
                "cp2k_input": inp,
                "input_status": "written",
            })
        except Exception as exc:
            updated.append({**rec, "input_status": f"input_gen_failed: {exc}"})
    return updated


def _bond_signature(struct: Structure) -> dict:
    max_cutoff = 5.2
    edges = []
    for i, site in enumerate(struct):
        ri = _get_covalent_radius(str(site.specie.symbol))
        for neigh in struct.get_neighbors(site, max_cutoff):
            j = int(neigh.index)
            if j <= i:
                continue
            rj = _get_covalent_radius(str(struct[j].specie.symbol))
            if neigh.nn_distance <= ri + rj + 0.5 and neigh.nn_distance > 0.35:
                a, b = sorted((str(site.specie.symbol), str(struct[j].specie.symbol)))
                edges.append((a, b))
    edge_counts = {}
    for key in edges:
        edge_counts[key] = edge_counts.get(key, 0) + 1
    degree_counts = {}
    degree = [0] * len(struct)
    for i, site in enumerate(struct):
        ri = _get_covalent_radius(str(site.specie.symbol))
        for neigh in struct.get_neighbors(site, max_cutoff):
            j = int(neigh.index)
            if j == i:
                continue
            rj = _get_covalent_radius(str(struct[j].specie.symbol))
            if neigh.nn_distance <= ri + rj + 0.5 and neigh.nn_distance > 0.35:
                degree[i] += 1
    for site, deg in zip(struct, degree):
        key = (str(site.specie.symbol), deg)
        degree_counts[key] = degree_counts.get(key, 0) + 1
    return {"edge_counts": edge_counts, "degree_counts": degree_counts}


def _rmsd_fractional_pbc(initial: Structure, relaxed: Structure) -> float | None:
    if len(initial) != len(relaxed):
        return None
    if [str(s.specie.symbol) for s in initial] != [str(s.specie.symbol) for s in relaxed]:
        return None
    dfrac = relaxed.frac_coords - initial.frac_coords
    dfrac -= np.round(dfrac)
    dcart = np.dot(dfrac, relaxed.lattice.matrix)
    return float(np.sqrt(np.mean(np.sum(dcart * dcart, axis=1))))


def _add_metrics(rec: dict, out_text: str, relaxed_cif: str | None) -> dict:
    energies = parse_energy_hartrees(out_text)
    initial_e = energies[0] if energies else None
    final_e = energies[-1] if energies else None
    updated = {
        **rec,
        "initial_e_hartree": initial_e,
        "final_e_hartree": final_e,
        "dft_e_hartree": final_e,
    }
    if initial_e is not None and final_e is not None:
        updated["delta_e_hartree"] = final_e - initial_e
        updated["delta_e_ev_per_atom"] = (
            (final_e - initial_e) * config.HARTREE_TO_EV / max(1, rec.get("n_atoms") or 1)
        )
    if relaxed_cif and Path(relaxed_cif).exists():
        initial = Structure.from_file(rec["clean_cif"])
        relaxed = Structure.from_file(relaxed_cif)
        updated["relaxed_cif"] = relaxed_cif
        updated["volume_change_pct"] = (
            (relaxed.lattice.volume - initial.lattice.volume) / initial.lattice.volume * 100.0
        )
        updated["rmsd_a"] = _rmsd_fractional_pbc(initial, relaxed)
        updated["bond_signature_retained"] = _bond_signature(initial) == _bond_signature(relaxed)
    return updated


def harvest(records: list[dict], results_dir: Path) -> list[dict]:
    updated = []
    for rec in records:
        run_dir = results_dir / rec["topology"] / "runs" / rec["id"]
        out_path = run_dir / "cp2k_run.out"
        relaxed_dir = results_dir / rec["topology"] / "relaxed"
        relaxed_dir.mkdir(parents=True, exist_ok=True)
        relaxed_cif = relaxed_dir / f"{rec['id']}.cif"
        if not out_path.exists():
            updated.append({**rec, "status": "no_output_file"})
            continue
        ok, detail = write_relaxed_cif_from_restart(str(run_dir), str(relaxed_cif))
        out_text = out_path.read_text(errors="replace")
        run_type = rec.get("relax_run_type") or config.RELAXATION_RUN_TYPE
        completion = parse_relaxation_completion(out_text, run_type)
        new_rec = {
            **rec,
            "relax_run_type": run_type,
            "relaxed_extract_status": detail,
            **(completion if ok else {}),
            "status": completion["status"] if ok else "relax_extract_failed",
        }
        new_rec = _add_metrics(new_rec, out_text, str(relaxed_cif) if ok else None)
        updated.append(new_rec)
    write_basin_pairs(updated, results_dir)
    return updated


def write_basin_pairs(records: list[dict], results_dir: Path) -> None:
    path = results_dir / "basin_pairs.csv"
    fields = ["id_a", "id_b", "same_formula", "same_n_atoms", "rmsd_a", "volume_delta_pct", "same_basin_proxy"]
    rows = []
    relaxed = [r for r in records if r.get("relaxed_cif") and Path(r["relaxed_cif"]).exists()]
    for i, a in enumerate(relaxed):
        sa = Structure.from_file(a["relaxed_cif"])
        for b in relaxed[i + 1:]:
            sb = Structure.from_file(b["relaxed_cif"])
            same_formula = sa.composition.reduced_formula == sb.composition.reduced_formula
            same_n_atoms = len(sa) == len(sb)
            rmsd = _rmsd_fractional_pbc(sa, sb) if same_n_atoms else None
            vol_delta = abs(sa.lattice.volume - sb.lattice.volume) / max(sa.lattice.volume, sb.lattice.volume) * 100.0
            rows.append({
                "id_a": a["id"],
                "id_b": b["id"],
                "same_formula": same_formula,
                "same_n_atoms": same_n_atoms,
                "rmsd_a": rmsd,
                "volume_delta_pct": vol_delta,
                "same_basin_proxy": bool(same_formula and same_n_atoms and rmsd is not None and rmsd < 0.25 and vol_delta < 2.0),
            })
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _run_records(records: list[dict], results_dir: Path, args) -> list[dict]:
    if args.profile == "local_safe":
        too_large = [r for r in records if (r.get("n_atoms") or 0) > args.local_max_atoms]
        if too_large:
            print("Local-safe profile will skip large structures:")
            for rec in too_large:
                print(f"  {rec['id']}: n_atoms={rec.get('n_atoms')}")
        records_to_run = [r for r in records if (r.get("n_atoms") or 0) <= args.local_max_atoms]
    else:
        records_to_run = list(records)

    if args.profile == "local_safe" and not args.yes_run_cp2k:
        print("Refusing to launch CP2K in local-safe mode without --yes-run-cp2k.")
        print("Use --write-inputs first, then run one ffc case only if the machine is idle.")
        return records

    if not Path(config.CP2K_BIN).exists():
        print(f"CP2K binary not found: {config.CP2K_BIN}")
        return [{**r, "status": "cp2k_bin_missing"} if r in records_to_run else r for r in records]

    old_omp = config.CP2K_OMP
    config.CP2K_OMP = config.RELAXATION_CP2K_OMP[args.profile]
    cutoff = args.cutoff_ry
    if args.profile == "local_safe" and cutoff is None:
        cutoff = config.SMOKE_CUTOFF_RY
    walltime = config.RELAXATION_WALLTIME_SECONDS[args.profile]
    max_mem = config.RELAXATION_MAX_MEM_GB[args.profile]
    try:
        by_id = {r["id"]: r for r in records}
        target_ids = {args.only} if args.only else {r["id"] for r in records_to_run}
        for rec in records_to_run:
            if rec["id"] not in target_ids:
                continue
            by_id[rec["id"]] = run_one(
                rec, str(results_dir), rec["topology"],
                relax_mode_override=args.run_type,
                walltime_seconds=walltime,
                cutoff_ry=cutoff,
                max_mem_gb=max_mem,
            )
        return [by_id[r["id"]] for r in records]
    finally:
        config.CP2K_OMP = old_omp


def main() -> None:
    parser = argparse.ArgumentParser(description="Small CP2K relaxation validation batch")
    parser.add_argument("--results-dir", default=config.RELAXATION_RESULTS_DIR)
    parser.add_argument("--preset", choices=["ffc-smoke", "fig2", "reviewer7", "reviewer9"], default="ffc-smoke")
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--write-inputs", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--harvest", action="store_true")
    parser.add_argument("--run-type", choices=["GEO_OPT", "CELL_OPT"], default=config.RELAXATION_RUN_TYPE)
    parser.add_argument("--profile", choices=["local_safe", "server"], default="local_safe")
    parser.add_argument("--cutoff-ry", type=int, default=None)
    parser.add_argument("--local-max-atoms", type=int, default=220)
    parser.add_argument("--only")
    parser.add_argument("--yes-run-cp2k", action="store_true")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if args.prepare or not _selection_path(results_dir).exists():
        records = prepare_records(args.preset, results_dir, only_id=args.only)
        save_records(records, results_dir)
        print(f"Prepared {len(records)} relaxation records in {results_dir}")
    else:
        records = load_records(results_dir)

    if args.write_inputs:
        cutoff = args.cutoff_ry
        if args.profile == "local_safe" and cutoff is None:
            cutoff = config.SMOKE_CUTOFF_RY
        records = write_inputs(records, results_dir, args.run_type, cutoff)
        save_records(records, results_dir)
        print(f"Wrote CP2K inputs for {len(records)} records")

    if args.run:
        records = _run_records(records, results_dir, args)
        save_records(records, results_dir)

    if args.harvest:
        records = harvest(records, results_dir)
        save_records(records, results_dir)
        print(f"Harvested relaxation outputs into {results_dir / 'relaxation_manifest.csv'}")


if __name__ == "__main__":
    main()
