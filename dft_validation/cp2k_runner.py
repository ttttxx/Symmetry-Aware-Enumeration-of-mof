import os, re, subprocess, json, shutil
from pymatgen.core import Structure
from dft_validation import config

OPEN_SHELL_METALS = {"Cu", "Ni", "Fe", "Co", "Mn", "V", "Cr", "Mo", "W"}
HARTREE_TO_EV = config.HARTREE_TO_EV


def _resolve_relax_mode(topology: str, n_atoms: int,
                         relax_mode_override: str | None = None) -> str:
    if relax_mode_override:
        return relax_mode_override
    rm = config.RELAX_MODE
    mode = rm[topology] if isinstance(rm, dict) else rm
    if mode != "auto":
        return mode
    if n_atoms < 100:
        return "CELL_OPT"
    if n_atoms <= 300:
        return "GEO_OPT"
    return "ENERGY"


def _is_open_shell(struct: Structure, multiplicity: int) -> bool:
    if multiplicity > 1:
        return True
    elems = {str(site.specie.symbol) for site in struct}
    return bool(elems & OPEN_SHELL_METALS)


def generate_cp2k_input(record: dict, run_dir: str, topology: str,
                         relax_mode_override: str | None = None) -> str:
    """Generate CP2K .inp file in run_dir. Returns path to .inp file."""
    clean_cif = record.get("clean_cif") or record.get("init_cif")
    struct = Structure.from_file(clean_cif)
    n_atoms = len(struct)
    multiplicity = record.get("multiplicity", 1)
    run_type = _resolve_relax_mode(topology, n_atoms, relax_mode_override)
    open_shell = _is_open_shell(struct, multiplicity)

    # Unique element list
    elements = sorted({str(site.specie.symbol) for site in struct})
    missing = [e for e in elements if e not in config.ELEMENT_KIND]
    if missing:
        raise ValueError(f"No basis/potential for elements: {missing}")

    # KIND blocks
    kind_blocks = ""
    for elem in elements:
        pot, basis = config.ELEMENT_KIND[elem]
        kind_blocks += f"""    &KIND {elem}
      BASIS_SET {basis}
      POTENTIAL {pot}
    &END KIND\n"""

    # SCF block
    if open_shell:
        scf_block = f"""      SCF_GUESS ATOMIC
      MAX_SCF 200
      EPS_SCF 1.0E-6
      ADDED_MOS 20
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE 300
      &END SMEAR
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.1
        NBROYDEN 8
      &END MIXING"""
    else:
        scf_block = """      SCF_GUESS ATOMIC
      MAX_SCF 200
      EPS_SCF 1.0E-6
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        MAX_SCF 10
        EPS_SCF 1.0E-5
      &END OUTER_SCF"""

    uks_line   = "      UKS .TRUE." if open_shell else ""
    multi_line = f"      MULTIPLICITY {multiplicity}" if multiplicity > 1 else ""

    # MOTION block (only for relaxation)
    motion_block = ""
    if run_type in ("GEO_OPT", "CELL_OPT"):
        motion_block = f"""&MOTION
  &{run_type}
    MAX_ITER 200
    MAX_FORCE 4.5E-4
    RMS_FORCE 3.0E-4
  &END {run_type}
&END MOTION
"""

    # Cartesian coords
    coord_lines = ""
    for site in struct:
        sym = str(site.specie.symbol)
        x, y, z = site.coords
        coord_lines += f"      {sym}  {x:.8f}  {y:.8f}  {z:.8f}\n"

    a, b, c = struct.lattice.abc
    al, be, ga = struct.lattice.angles

    inp = f"""&GLOBAL
  PROJECT cp2k_run
  RUN_TYPE {run_type}
  PRINT_LEVEL LOW
&END GLOBAL

{motion_block}&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME {config.CP2K_DATA_DIR}/{config.BASIS_SET_FILE}
    POTENTIAL_FILE_NAME {config.CP2K_DATA_DIR}/GTH_POTENTIALS
    CHARGE {config.CHARGE}
{uks_line}
{multi_line}
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
    &MGRID
      CUTOFF {config.PW_CUTOFF_RY}
      REL_CUTOFF {config.REL_CUTOFF_RY}
    &END MGRID
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          PARAMETER_FILE_NAME {config.CP2K_DATA_DIR}/dftd3.dat
          REFERENCE_FUNCTIONAL PBE
          R_CUTOFF 15.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
    &SCF
{scf_block}
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      A {a:.6f} 0.0 0.0
      B 0.0 {b:.6f} 0.0
      C 0.0 0.0 {c:.6f}
      ALPHA_BETA_GAMMA {al:.4f} {be:.4f} {ga:.4f}
      PERIODIC XYZ
    &END CELL
    &COORD
      SCALED .FALSE.
{coord_lines}    &END COORD
{kind_blocks}  &END SUBSYS
&END FORCE_EVAL
"""

    inp_path = os.path.join(run_dir, "cp2k_run.inp")
    with open(inp_path, "w") as f:
        f.write(inp)
    return inp_path


def parse_energy(output_text: str, n_atoms: int) -> float | None:
    """
    Extract last 'ENERGY| Total FORCE_EVAL' line, convert Hartree → eV/atom.
    Returns None if not found.
    """
    matches = re.findall(
        r"ENERGY\| Total FORCE_EVAL[^\n]*\s+(-?\d+\.\d+)", output_text
    )
    if not matches:
        return None
    e_hartree = float(matches[-1])
    return e_hartree * HARTREE_TO_EV / n_atoms


def run_one(record: dict, results_dir: str, topology: str,
            relax_mode_override: str | None = None,
            walltime_seconds: int = 3600) -> dict:
    """
    Run CP2K for one record. Returns updated record with dft_e_per_atom and status.
    Skips if relaxed CIF already exists (resumable).
    """
    rec_id = record["id"]
    run_dir = os.path.join(results_dir, topology, "runs", rec_id)
    relaxed_dir = os.path.join(results_dir, topology, "relaxed")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(relaxed_dir, exist_ok=True)

    relaxed_cif = os.path.join(relaxed_dir, f"{rec_id}.cif")
    if os.path.exists(relaxed_cif):
        print(f"  [runner] {rec_id}: already done, skipping")
        return {**record, "status": "done_cached"}

    if record.get("status") not in ("clean", "preclean_done", "selected", "assembled"):
        print(f"  [runner] {rec_id}: skipping (status={record.get('status')})")
        return record

    try:
        inp_path = generate_cp2k_input(record, run_dir, topology, relax_mode_override)
    except Exception as e:
        return {**record, "status": f"input_gen_failed: {e}"}

    out_path = os.path.join(run_dir, "cp2k_run.out")
    abs_inp = os.path.abspath(inp_path)
    abs_out = os.path.abspath(out_path)
    cmd = f"{config.CP2K_RUN_CMD} -i {abs_inp} -o {abs_out}"
    print(f"  [runner] {rec_id}: running — {cmd}")

    try:
        result = subprocess.run(
            cmd, shell=True, timeout=walltime_seconds,
            capture_output=True, text=True, cwd=run_dir
        )
    except subprocess.TimeoutExpired:
        return {**record, "status": "timeout"}
    except Exception as e:
        return {**record, "status": f"subprocess_error: {e}"}

    if not os.path.exists(out_path):
        return {**record, "status": "no_output_file"}

    out_text = open(out_path).read()
    n_atoms = record.get("n_atoms") or 1
    e_per_atom = parse_energy(out_text, n_atoms)

    if e_per_atom is None:
        snippet = out_text[-500:]
        return {**record, "status": "parse_failed", "cp2k_log_tail": snippet}

    run_type = _resolve_relax_mode(topology, n_atoms, relax_mode_override)
    src = record.get("clean_cif") or record.get("init_cif")
    if src:
        shutil.copy(src, relaxed_cif)

    e_hartree_matches = re.findall(
        r"ENERGY\| Total FORCE_EVAL[^\n]*\s+(-?\d+\.\d+)", out_text
    )
    e_hartree = float(e_hartree_matches[-1]) if e_hartree_matches else None

    return {**record,
            "dft_e_per_atom": e_per_atom,
            "dft_e_hartree": e_hartree,
            "relaxed_cif": relaxed_cif,
            "status": "converged"}
