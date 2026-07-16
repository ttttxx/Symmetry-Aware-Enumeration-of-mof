import os, json, sys
import numpy as np
from pymatgen.core import Structure
from dft_validation import config

def check_element_coverage(struct: Structure) -> list:
    """Return list of element symbols not in config.ELEMENT_KIND."""
    known = set(config.ELEMENT_KIND.keys())
    present = {str(sp.symbol) for sp in struct.composition.elements}
    return sorted(present - known)

def compute_electron_parity(struct: Structure) -> tuple:
    """
    Sum valence electrons (from config.ELEMENT_VALENCE) minus config.CHARGE.
    Returns (n_valence_electrons, multiplicity).
    multiplicity = 2 if n_electrons is odd, else 1.
    """
    n_e = sum(
        config.ELEMENT_VALENCE.get(str(site.specie.symbol), 0)
        for site in struct
    ) - config.CHARGE
    mult = 2 if n_e % 2 == 1 else 1
    return n_e, mult

def min_nonbond_distance(struct: Structure) -> float:
    """Return the minimum pairwise distance between all sites."""
    min_d = float("inf")
    for i in range(len(struct)):
        for j in range(i + 1, len(struct)):
            d = struct.get_distance(i, j)
            if d < min_d:
                min_d = d
    return min_d

def preclean_xtb(struct: Structure) -> Structure | None:
    """
    Run GFN1-xTB FIRE geometry optimisation (periodic).
    Returns optimised Structure or None on failure.
    Energy is discarded — geometry only.
    """
    try:
        from xtb.ase.calculator import XTB
        from ase.optimize import FIRE
        from pymatgen.io.ase import AseAtomsAdaptor

        atoms = AseAtomsAdaptor.get_atoms(struct)
        atoms.calc = XTB(method="GFN1-xTB")
        opt = FIRE(atoms, logfile=None)
        opt.run(fmax=config.PRECLEAN_FMAX, steps=config.PRECLEAN_STEPS)
        return AseAtomsAdaptor.get_structure(atoms)
    except Exception as e:
        print(f"  [sanitize] xTB failed: {e}")
        return None

def _preclean(struct: Structure) -> Structure | None:
    """Dispatch to the configured pre-clean engine."""
    engine = config.PRECLEAN_ENGINE
    if engine == "xtb":
        return preclean_xtb(struct)
    if engine == "chgnet":
        try:
            from chgnet.model import StructOptimizer
            result = StructOptimizer().relax(struct, relax_cell=False)
            return result["final_structure"]
        except Exception as e:
            print(f"  [sanitize] CHGNet failed: {e}")
            return None
    return struct  # "none" engine — return as-is

def sanitize_one(record: dict, results_dir: str) -> dict:
    """
    Run the full sanitize pipeline for a single record.
    Updates record fields: status, clean_cif, n_electrons, multiplicity.
    """
    if record.get("status") not in ("assembled", "selected"):
        return record  # skip already-failed or already-clean

    init_cif = record.get("init_cif")
    if not init_cif or not os.path.exists(init_cif):
        return {**record, "status": "missing_init_cif"}

    topology = record["topology"]
    clean_dir = os.path.join(results_dir, topology, "clean")
    os.makedirs(clean_dir, exist_ok=True)

    try:
        struct = Structure.from_file(init_cif)
    except Exception as e:
        return {**record, "status": f"load_failed: {e}"}

    # 1. Wrap fractional coords
    struct = Structure(struct.lattice, struct.species,
                       struct.frac_coords % 1.0,
                       coords_are_cartesian=False)

    # 2. Element coverage check
    missing = check_element_coverage(struct)
    if missing:
        print(f"  [sanitize] {record['id']}: missing kinds {missing}")
        return {**record, "status": f"missing_kind:{','.join(missing)}"}

    # 3. Electron parity
    n_e, mult = compute_electron_parity(struct)
    if mult == 2:
        print(f"  [sanitize] {record['id']}: odd electron count ({n_e}) → doublet")

    # 4. Min-distance check → pre-relax if needed
    min_d = min_nonbond_distance(struct)
    status = "clean"
    if min_d < config.MIN_DIST_FLOOR:
        print(f"  [sanitize] {record['id']}: min_dist={min_d:.3f} Å < {config.MIN_DIST_FLOOR} — pre-relaxing")
        relaxed = _preclean(struct)
        if relaxed is None:
            return {**record, "status": "preclean_failed",
                    "n_electrons": n_e, "multiplicity": mult}
        min_d2 = min_nonbond_distance(relaxed)
        if min_d2 < config.CLASH_HARD_LIMIT:
            print(f"  [sanitize] {record['id']}: still clashing ({min_d2:.3f} Å) — discarding")
            return {**record, "status": "clash_unfixable",
                    "n_electrons": n_e, "multiplicity": mult}
        struct = relaxed
        status = "preclean_done"

    # Save clean CIF
    cif_name = os.path.basename(init_cif)
    clean_path = os.path.join(clean_dir, cif_name)
    struct.to(filename=clean_path)

    return {**record,
            "status": status,
            "clean_cif": clean_path,
            "n_electrons": n_e,
            "multiplicity": mult,
            "n_atoms": len(struct)}

def sanitize_all(records: list, results_dir: str) -> list:
    """Sanitize every record; write updated JSON; return updated list."""
    updated = []
    for rec in records:
        print(f"  [sanitize] Processing {rec['id']} …")
        updated.append(sanitize_one(rec, results_dir))

    json_path = os.path.join(results_dir, "selected_structures.json")
    with open(json_path, "w") as f:
        json.dump(updated, f, indent=2)

    ok = sum(1 for r in updated if r.get("status") in ("clean", "preclean_done"))
    print(f"[sanitize] Done: {ok}/{len(updated)} structures clean/pre-cleaned")
    return updated
