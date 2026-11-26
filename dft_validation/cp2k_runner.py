import os, re, subprocess, json, shutil, time, threading, signal, glob
import psutil
from pymatgen.core import Lattice, Structure
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

def _run_with_oom_guard(cmd_list: list, env: dict, cwd: str,
                         timeout_sec: int, max_mem_gb: float) -> str | None:
    """
    Run cmd_list as a new process group; kill the entire group if the
    cumulative RSS of parent + all recursive children exceeds max_mem_gb.
    Polls every 2 seconds. Returns None on clean exit, or a status string
    ('failed_oom', 'timeout', 'popen_error: ...') on failure.
    """
    kill_reason: list[str | None] = [None]
    kill_lock = threading.Lock()

    try:
        proc = subprocess.Popen(cmd_list, env=env, cwd=cwd,
                                preexec_fn=os.setsid)
    except Exception as exc:
        return f"popen_error: {exc}"

    def _kill_process_group(reason: str, sig: int = signal.SIGTERM) -> None:
        with kill_lock:
            if kill_reason[0] is None:
                kill_reason[0] = reason
            if proc.poll() is not None:
                return
            try:
                pgid = os.getpgid(proc.pid)
                os.killpg(pgid, sig)
            except ProcessLookupError:
                return

    def _monitor():
        while proc.poll() is None:
            try:
                parent = psutil.Process(proc.pid)
                procs = [parent] + parent.children(recursive=True)
                total_rss = sum(p.memory_info().rss for p in procs
                                if p.is_running())
                if total_rss / (1024 ** 3) > max_mem_gb:
                    _kill_process_group("failed_oom", signal.SIGKILL)
                    return
            except (psutil.NoSuchProcess, ProcessLookupError):
                return
            time.sleep(2)

    monitor = threading.Thread(target=_monitor, daemon=True)
    monitor.start()

    old_handlers = {}

    def _handle_stop(signum, _frame):
        _kill_process_group(f"interrupted_signal_{signum}", signal.SIGTERM)
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            _kill_process_group(f"interrupted_signal_{signum}", signal.SIGKILL)
            proc.wait()
        raise KeyboardInterrupt(f"received signal {signum}; stopped CP2K")

    if threading.current_thread() is threading.main_thread():
        for sig in (signal.SIGINT, signal.SIGTERM):
            old_handlers[sig] = signal.getsignal(sig)
            signal.signal(sig, _handle_stop)

    try:
        proc.wait(timeout=timeout_sec)
    except subprocess.TimeoutExpired:
        _kill_process_group("timeout", signal.SIGKILL)
        proc.wait()
    except BaseException:
        _kill_process_group("interrupted", signal.SIGTERM)
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            _kill_process_group("interrupted", signal.SIGKILL)
            proc.wait()
        raise
    finally:
        for sig, old_handler in old_handlers.items():
            signal.signal(sig, old_handler)

    monitor.join(timeout=5)
    return kill_reason[0]

def generate_cp2k_input(record: dict, run_dir: str, topology: str,
                         relax_mode_override: str | None = None,
                         cutoff_ry: int | None = None) -> str:
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

    # SCF WFN restart: CP2K auto-finds cp2k_run-RESTART.wfn in cwd.
    restart_wfn = os.path.join(run_dir, "cp2k_run-RESTART.wfn")
    wfn_exists  = os.path.exists(restart_wfn)
    scf_guess   = "RESTART" if wfn_exists else "ATOMIC"
    if wfn_exists:
        print(f"  [cp2k_inp] WFN restart detected — using SCF_GUESS RESTART")

    # GEO_OPT / CELL_OPT coordinate restart: use the newest numeric restart.
    geo_restart = _latest_restart_file(run_dir, include_backups=False)
    if geo_restart:
        print(f"  [cp2k_inp] GEO restart detected — {os.path.basename(geo_restart)}")
    ext_restart_block = (
        f"&EXT_RESTART\n  RESTART_FILE_NAME {os.path.abspath(geo_restart)}\n&END EXT_RESTART\n\n"
        if geo_restart else ""
    )

    # SCF block
    if open_shell:
        scf_block = f"""      SCF_GUESS {scf_guess}
      MAX_SCF 300
      EPS_SCF 5.0E-5
      ADDED_MOS 30
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE 300
      &END SMEAR
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.1
        NBROYDEN 8
      &END MIXING
      IGNORE_CONVERGENCE_FAILURE .TRUE."""
    else:
        scf_block = f"""      SCF_GUESS {scf_guess}
      MAX_SCF 300
      EPS_SCF 5.0E-5
      ADDED_MOS 30
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE 300
      &END SMEAR
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.2
        NBROYDEN 8
      &END MIXING
      IGNORE_CONVERGENCE_FAILURE .TRUE."""

    uks_line   = "      UKS .TRUE." if open_shell else ""
    multi_line = f"      MULTIPLICITY {multiplicity}" if multiplicity > 1 else ""

    # MOTION block (only for relaxation)
    motion_block = ""
    if run_type in ("GEO_OPT", "CELL_OPT"):
        motion_block = f"""&MOTION
  &{run_type}
    MAX_ITER 200
    MAX_FORCE 1.0E-3
    RMS_FORCE 7.0E-4
  &END {run_type}
  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      ADD_LAST NUMERIC
      &EACH
        GEO_OPT 1
        CELL_OPT 1
      &END EACH
    &END TRAJECTORY
    &CELL
      &EACH
        GEO_OPT 1
        CELL_OPT 1
      &END EACH
    &END CELL
    &RESTART
      BACKUP_COPIES 1
      &EACH
        GEO_OPT 1
        CELL_OPT 1
      &END EACH
    &END RESTART
  &END PRINT
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

    pw_cut = cutoff_ry if cutoff_ry is not None else config.PW_CUTOFF_RY
    rel_cut = max(30, int(config.REL_CUTOFF_RY * pw_cut / config.PW_CUTOFF_RY))
    stress_line = "  STRESS_TENSOR ANALYTICAL\n" if run_type == "CELL_OPT" else ""

    inp = f"""{ext_restart_block}&GLOBAL
  PROJECT cp2k_run
  RUN_TYPE {run_type}
  PRINT_LEVEL LOW
&END GLOBAL

{motion_block}&FORCE_EVAL
  METHOD Quickstep
{stress_line}
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
      CUTOFF {pw_cut}
      REL_CUTOFF {rel_cut}
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

def parse_energy_hartrees(output_text: str) -> list[float]:
    """Return all CP2K total FORCE_EVAL energies in Hartree."""
    return [
        float(x)
        for x in re.findall(
            r"ENERGY\| Total FORCE_EVAL[^\n]*\s+(-?\d+\.\d+)",
            output_text,
        )
    ]

def parse_relaxation_completion(output_text: str, run_type: str | None = None) -> dict:
    """Classify CP2K relaxation output without treating restart harvest as convergence."""
    upper = output_text.upper()
    completed_markers = [
        "OPTIMIZATION COMPLETED",
        "GEOMETRY OPTIMIZATION COMPLETED",
        "CELL OPTIMIZATION COMPLETED",
    ]
    program_ended = "PROGRAM ENDED" in upper
    optimization_completed = any(marker in upper for marker in completed_markers)
    steps = [
        int(x)
        for x in re.findall(r"OPTIMIZATION STEP:\s*(\d+)", output_text, flags=re.IGNORECASE)
    ]
    criteria = re.findall(
        r"OPT\|\s*([^|\n]*?)\s+is converged\s+(YES|NO)",
        output_text,
        flags=re.IGNORECASE,
    )
    last_flags = {}
    for label, flag in criteria:
        label = re.sub(r"\s+", " ", label.strip()).lower().replace(" ", "_")
        last_flags[label] = flag.upper()
    all_flags_yes = None
    if last_flags:
        all_flags_yes = all(flag == "YES" for flag in last_flags.values())

    status = "relaxation_converged"
    if run_type in {"GEO_OPT", "CELL_OPT"}:
        status = f"{run_type.lower()}_converged"
    if not (program_ended and optimization_completed):
        status = "harvested_partial_restart"

    return {
        "status": status,
        "cp2k_program_ended": program_ended,
        "relaxation_completed_marker": optimization_completed,
        "optimization_steps": steps[-1] if steps else None,
        "relaxation_convergence_flags_all_yes": all_flags_yes,
        "relaxation_convergence_flags": ";".join(
            f"{label}={flag}" for label, flag in sorted(last_flags.items())
        ),
    }

def _numbers_from_cp2k_line(parts: list[str]) -> list[float]:
    nums = []
    for token in parts:
        clean = token.replace("D", "E")
        try:
            nums.append(float(clean))
        except ValueError:
            continue
    return nums

def _restart_sort_key(path: str) -> tuple[int, int, float]:
    name = os.path.basename(path)
    match = re.search(r"_(\d+)\.restart(?:\.bak-\d+)?$", name)
    step = int(match.group(1)) if match else -1
    backup_penalty = 0 if ".bak-" not in name else -1
    try:
        mtime = os.path.getmtime(path)
    except OSError:
        mtime = 0.0
    return (step, backup_penalty, mtime)

def _latest_restart_file(run_dir: str, include_backups: bool = True) -> str | None:
    candidates = []
    patterns = ["*.restart"]
    if include_backups:
        patterns.append("*.restart.bak-*")
    for pattern in patterns:
        candidates.extend(glob.glob(os.path.join(run_dir, pattern)))
    if not candidates:
        return None
    return max(candidates, key=_restart_sort_key)

def _parse_structure_from_restart(restart_path: str) -> Structure:
    """Parse final cell and coordinates from a CP2K restart-style input file."""
    lines = open(restart_path, errors="replace").read().splitlines()
    cell = {}
    species = []
    coords = []
    scaled = False
    in_cell = False
    in_coord = False

    for raw in lines:
        line = raw.strip()
        upper = line.upper()
        if not line or upper.startswith(("!", "#", "@")):
            continue
        if upper.startswith("&CELL"):
            in_cell = True
            continue
        if in_cell and upper.startswith("&END"):
            in_cell = False
            continue
        if upper.startswith("&COORD"):
            in_coord = True
            continue
        if in_coord and upper.startswith("&END"):
            in_coord = False
            continue

        if in_cell:
            parts = line.split()
            key = parts[0].upper()
            if key in {"A", "B", "C"}:
                nums = _numbers_from_cp2k_line(parts[1:])
                if len(nums) >= 3:
                    cell[key] = nums[:3]
            continue

        if in_coord:
            parts = line.split()
            key = parts[0].upper()
            if key == "SCALED":
                if len(parts) > 1:
                    scaled = parts[1].strip().upper().startswith(".T")
                continue
            if key == "UNIT":
                continue
            nums = _numbers_from_cp2k_line(parts[1:])
            if len(nums) >= 3:
                species.append(parts[0])
                coords.append(nums[:3])

    missing = [k for k in ("A", "B", "C") if k not in cell]
    if missing:
        raise ValueError(f"restart file lacks cell vectors: {missing}")
    if not species:
        raise ValueError("restart file lacks coordinates")

    lattice = Lattice([cell["A"], cell["B"], cell["C"]])
    return Structure(lattice, species, coords, coords_are_cartesian=not scaled)

def write_relaxed_cif_from_restart(run_dir: str, relaxed_cif: str) -> tuple[bool, str]:
    """Write a real relaxed CIF from the latest CP2K restart file."""
    restart = _latest_restart_file(run_dir)
    if restart is None:
        return False, "no_restart_file"
    try:
        struct = _parse_structure_from_restart(restart)
        struct.to(filename=relaxed_cif)
    except Exception as exc:
        return False, f"restart_parse_failed: {exc}"
    return True, os.path.basename(restart)

def run_one(record: dict, results_dir: str, topology: str,
            relax_mode_override: str | None = None,
            walltime_seconds: int = 3600,
            cutoff_ry: int | None = None,
            max_mem_gb: float | None = None) -> dict:
    """
    Run CP2K for one record. Returns updated record with dft_e_per_atom and status.
    Skips if relaxed CIF already exists (resumable).
    Uses 1 MPI rank + CP2K_OMP OpenMP threads; OMP_NUM_THREADS injected via env.
    Process tree RSS is monitored every 2 s; exceeding max_mem_gb kills the group.
    """
    rec_id = record["id"]
    run_dir = os.path.join(results_dir, topology, "runs", rec_id)
    relaxed_dir = os.path.join(results_dir, topology, "relaxed")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(relaxed_dir, exist_ok=True)

    relaxed_cif = os.path.join(relaxed_dir, f"{rec_id}.cif")
    if os.path.exists(relaxed_cif):
        run_type = _resolve_relax_mode(
            topology, record.get("n_atoms") or 1, relax_mode_override
        )
        out_path = os.path.join(run_dir, "cp2k_run.out")
        if run_type not in ("GEO_OPT", "CELL_OPT"):
            print(f"  [runner] {rec_id}: already done, skipping")
            return {**record, "status": "done_cached"}
        if os.path.exists(out_path):
            out_text = open(out_path, errors="replace").read()
            completion = parse_relaxation_completion(out_text, run_type)
            if completion["status"].endswith("_converged"):
                print(f"  [runner] {rec_id}: converged relaxation cached, skipping")
                return {**record, **completion, "status": "done_cached"}
        print(f"  [runner] {rec_id}: partial relaxed CIF exists; attempting restart")

    if record.get("status") not in (
        "clean", "preclean_done", "selected", "assembled",
        "cp2k_bin_missing", "no_output_file", "failed_oom", "timeout",
        "harvested_partial_restart", "relax_extract_failed", "parse_failed",
    ):
        print(f"  [runner] {rec_id}: skipping (status={record.get('status')})")
        return record

    try:
        inp_path = generate_cp2k_input(record, run_dir, topology,
                                       relax_mode_override, cutoff_ry)
    except Exception as e:
        return {**record, "status": f"input_gen_failed: {e}"}

    out_path = os.path.join(run_dir, "cp2k_run.out")
    abs_inp  = os.path.abspath(inp_path)
    abs_out  = os.path.abspath(out_path)

    # Build env: inherit PATH etc., inject OMP threads explicitly.
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(config.CP2K_OMP)

    # Single-node static binary (.ssmp): OpenMP only, no MPI launcher.
    cmd_list = [config.CP2K_BIN, "-i", abs_inp, "-o", abs_out]
    print(f"  [runner] {rec_id}: running — {' '.join(cmd_list)}")
    print(f"  [runner] {rec_id}: output  — {abs_out}")

    fail_reason = _run_with_oom_guard(
        cmd_list, env=env, cwd=run_dir,
        timeout_sec=walltime_seconds,
        max_mem_gb=max_mem_gb if max_mem_gb is not None else config.MAX_MEM_GB,
    )
    if fail_reason:
        return {**record, "status": fail_reason}

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
    relaxed_extract_status = None
    if run_type in ("GEO_OPT", "CELL_OPT"):
        ok, detail = write_relaxed_cif_from_restart(run_dir, relaxed_cif)
        relaxed_extract_status = detail
        if not ok:
            e_hartrees = parse_energy_hartrees(out_text)
            e_hartree = e_hartrees[-1] if e_hartrees else None
            return {**record,
                    "dft_e_per_atom": e_per_atom,
                    "dft_e_hartree": e_hartree,
                    "relax_run_type": run_type,
                    "relaxed_extract_status": relaxed_extract_status,
                    "status": "relax_extract_failed"}
    elif src:
        shutil.copy(src, relaxed_cif)

    e_hartrees = parse_energy_hartrees(out_text)
    e_hartree = e_hartrees[-1] if e_hartrees else None
    completion = {}
    status = "converged"
    if run_type in ("GEO_OPT", "CELL_OPT"):
        completion = parse_relaxation_completion(out_text, run_type)
        status = completion["status"]

    return {**record,
            "dft_e_per_atom": e_per_atom,
            "dft_e_hartree": e_hartree,
            "relax_run_type": run_type,
            "relaxed_extract_status": relaxed_extract_status,
            "relaxed_cif": relaxed_cif,
            **completion,
            "status": status}
