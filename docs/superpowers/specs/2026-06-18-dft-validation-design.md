# DFT Validation Pipeline — Design Spec

**Date:** 2026-06-18  
**Goal:** Enumerate MOF structures by symmetry class (High/Med/Low), sanitize geometries, run CP2K PBE-D3(BJ) single-point or relaxation, and test whether DFT energy correlates inversely with symmetry index ρ (Spearman ρ_s per topology).

---

## 1. Environment (verified)

| Item | Value |
|------|-------|
| CP2K binary | `/opt/homebrew/bin/cp2k.psmp` (v2026.1, MPI+OMP) |
| CP2K_RUN_CMD | `mpirun -np 4 cp2k.psmp` |
| CP2K_DATA_DIR | `/opt/homebrew/Cellar/cp2k/2026.1_1/share/cp2k/data` |
| GFN1-xTB via CP2K | **Not compiled** (no `xtb` flag) |
| Standalone xTB | `xtb-python` + ASE at `/opt/anaconda3` (confirmed OK) |

---

## 2. Local vs Server Profiles (§7 decisions)

### A. Local smoke-test profile

**Primary run:** `USE_CLOSED_SHELL_VARIANT=True` — swap all `4c_Cu_1_Ch` → `4c_Cd_1_Ch`.  
Cd is closed-shell, uses OT SCF, converges reliably on 16 GB RAM. Validates the full pipeline.

**Cu probe:** Additionally select 1 High-symmetry ffc structure with Cu retained (ignoring the Cd override for this probe only). Run as `ENERGY` single-point to exercise the open-shell code path:
- `SPIN_POLARIZED=True`, `SMEARING="auto"` → Fermi-Dirac 300 K, diagonal SCF, UKS
- Lives alongside the Cd batch; a Cu SCF failure sets `status=failed` without stopping the batch.
- Config flag: `CU_PROBE = True`, `CU_PROBE_N = 1` (default on locally, default off on server since the full Cu run is the server's job).

### B. fjh strategy

**Run type:** `RELAX_MODE` is per-topology: `{"ffc": "auto", "fjh": "ENERGY"}`.  
fjh has ~800 atoms — single-point only, even on a server.

**Sample size:** `N_PER_CLASS` applies to both, but fjh is gated:

**Gating logic:** `FJH_GATE_ON_FFC = True` (default). After ffc completes, check:
- Spearman p < 0.1 (configurable: `FJH_GATE_SPEARMAN_P = 0.1`), AND
- Lowest-energy structure is in the High class.

If both pass → proceed with fjh. If either fails → print a warning and skip fjh (not an error). Gate can be disabled with `FJH_GATE_ON_FFC = False` for forced server runs.

---

## 3. Pre-relaxation Engine

CP2K's built-in xTB is unavailable. `PRECLEAN_ENGINE = "xtb"` routes to **standalone xtb-python** (ASE interface), not CP2K:

```python
from xtb.ase.calculator import XTB
from ase.optimize import BFGS
from ase.io import read, write

calc = XTB(method="GFN1-xTB")
atoms.calc = calc
opt = BFGS(atoms, logfile=None)
opt.run(fmax=0.1)   # loose convergence — geometry only
```

This is consistent across ffc and fjh (same engine regardless of atom count).  
Fallback chain: `"xtb"` → `"chgnet"` → `"none"`.

---

## 4. Package Structure

```
dft_validation/
├── __init__.py
├── config.py                # all tunable knobs
├── select_structures.py     # enumerate → bucket → sample → assemble CIFs
├── sanitize.py              # wrap → element check → clash → pre-relax
├── cp2k_runner.py           # generate input → run → parse energy
├── analyze.py               # post-relax symmetry + Spearman + figures
└── run_dft_validation.py    # orchestrator: select→sanitize→run→analyze
```

No existing business files are modified except the two compat fixes in §1 of task.txt, handled in-package.

---

## 5. config.py (final values)

```python
# --- Targets ---
USE_CLOSED_SHELL_VARIANT = True   # True → swap 4c_Cu_1_Ch → 4c_Cd_1_Ch
CU_PROBE      = True              # local only: run 1 Cu High-sym structure as open-shell test
CU_PROBE_N    = 1

VALIDATION_TARGETS = [
    {"topology": "ffc",
     "atom_choices": {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]},
     "forced_edge": "1B_4H_Ch"},
    {"topology": "fjh",
     "atom_choices": {3: ["3c_B_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]},
     "forced_edge": "1B_4H_Ch"},
]

# --- Sampling ---
N_PER_CLASS  = 2
SYM_CLASSES  = ["High", "Med", "Low"]
RANDOM_SEED  = 42
MAX_ATOMS    = 250

# --- CHGNet (default off) ---
USE_CHGNET      = False
CHGNET_PRECLEAN = False   # used only if PRECLEAN_ENGINE=="chgnet"
CHGNET_COMPARE  = False

# --- CP2K ---
CP2K_RUN_CMD  = "mpirun -np 4 cp2k.psmp"
CP2K_DATA_DIR = "/opt/homebrew/Cellar/cp2k/2026.1_1/share/cp2k/data"

# --- Electronic structure ---
SPIN_POLARIZED = "auto"
DFT_PLUS_U     = None
SMEARING       = "auto"

# --- DFT spec ---
XC_FUNCTIONAL  = "PBE"
DISPERSION     = "DFTD3(BJ)"
BASIS_SET      = "DZVP-MOLOPT-SR-GTH"
PW_CUTOFF_RY   = 600
REL_CUTOFF_RY  = 60
RELAX_MODE     = {"ffc": "auto", "fjh": "ENERGY"}   # per-topology dict or single string

# --- fjh gate ---
FJH_GATE_ON_FFC       = True
FJH_GATE_SPEARMAN_P   = 0.1

# --- Sanitize ---
MIN_DIST_FLOOR  = 0.9
PRECLEAN_ENGINE = "xtb"   # "xtb" → xtb-python+ASE | "chgnet" | "none"
SYMPREC_POSTRELAX = 0.3

RESULTS_DIR = "results/dft_validation"
```

---

## 6. Module Responsibilities

### select_structures.py
- Enumerates per §1 of task.txt (enrich → topology_requirements → MOFEnumerator with degree constraints).
- Buckets results using the exact logic from `physical_validation.py` (lines 105–120): top-2 stab values → High, interior >1 values → Med, stab==1 → Low.
- Samples `N_PER_CLASS` per bucket with `RANDOM_SEED`; skips any structure with `n_atoms > MAX_ATOMS`.
- If `USE_CLOSED_SHELL_VARIANT`: replaces `4c_Cu_1_Ch` → `4c_Cd_1_Ch` in `atom_choices` before enumeration.
- Cu probe: after the Cd batch, separately pick the 1 highest-ρ ffc result with original Cu atom_choices.
- Assembles CIFs via `StructureAssembler.assemble_and_save`; saves `selected_structures.json`.

### sanitize.py
- For each init structure, in order:
  1. pymatgen: wrap fractional coords, verify 3D periodic, non-degenerate cell.
  2. Element coverage check: every element in the structure must have a `DZVP-MOLOPT-SR-GTH` entry in `BASIS_MOLOPT` and a `GTH-PBE` entry in `GTH_POTENTIALS`; missing → `status=missing_kind`, skip.
  3. Minimum-distance check: if `min_dist < MIN_DIST_FLOOR` → pre-relax. After pre-relax, if `min_dist < 0.6 Å` → `status=clash_unfixable`, discard.
- `PRECLEAN_ENGINE="xtb"`: uses `xtb-python` ASE calculator with `GFN1-xTB`, `fmax=0.1` (geometry only, energy discarded).
- Saves clean CIFs; updates `selected_structures.json` with `status` and `clean_cif` path.

### cp2k_runner.py
- Generates CP2K input from clean geometry. Key choices:
  - `SPIN_POLARIZED="auto"`: auto-detects Cu/Ni/Fe/Co/Mn/V → `UKS True`, estimates `MULTIPLICITY` from unpaired metal count; otherwise closed-shell.
  - `SMEARING="auto"`: closed-shell insulator → OT (`DIIS` + `FULL_SINGLE_INVERSE`) + OUTER_SCF; open-shell/metallic → diagonal + Fermi-Dirac 300 K + `ADDED_MOS 20` + Broyden mixing.
  - `RELAX_MODE`: resolved from config (per-topology dict or string); `auto` maps `<100 atoms→CELL_OPT`, `100–300→GEO_OPT`, `>300→ENERGY`.
  - Γ-point only for large cells.
  - `SCF_GUESS ATOMIC`, never RESTART on first run.
- Runs `CP2K_RUN_CMD` in a per-structure subdirectory; wraps in try/except; failed structures get `status=failed` + log entry; batch continues.
- Parses last `ENERGY| Total FORCE_EVAL` line → `dft_e_per_atom = E / n_atoms`; exports relaxed structure.
- Resumable: skips structures where `relaxed/*.cif` already exists.

### analyze.py
- Re-runs `SpacegroupAnalyzer(symprec=SYMPREC_POSTRELAX)` on relaxed structures.
- Computes `scipy.stats.spearmanr(rho_values, dft_e_per_atom_values)` per topology (converged structures only).
- Writes `dft_vs_symmetry.csv`.
- Produces ρ–E scatter plot per topology (sym_class colored), dpi=300.
- If `CHGNET_COMPARE`: adds CHGNet-vs-DFT scatter with diagonal + ρ_s annotation.
- **fjh gate evaluation**: after ffc analyze, returns `(rho_s, p_val, lowest_in_high)` for the orchestrator to check.

### run_dft_validation.py
- `main()`: runs select→sanitize→run→analyze for ffc.
- If `FJH_GATE_ON_FFC`: evaluates gate; prints pass/fail + reason; skips fjh on fail.
- If fjh proceeds: same pipeline with topology-specific overrides.
- Terminal summary line: topology, Spearman ρ_s, p-value, n_converged/N_total, n_skipped/n_discarded, output paths.

---

## 7. Compatibility Fixes (in-package, minimal)

- **NetworkX v3 `shortest_path` returns generator**: in `sanitize.py` or wherever `periodic_structure_analysis` is imported, apply `import networkx as nx; _orig = nx.shortest_path; nx.shortest_path = lambda *a, **kw: dict(_orig(*a, **kw))` monkeypatch before first import of `external_lib`.
- **`enrich_graph_data_with_topology` must be called before assembly**: enforced by select_structures.py calling it unconditionally (follows main.py pattern).

---

## 8. Acceptance

- `python dft_validation/run_dft_validation.py` runs end-to-end on local profile (ffc, N_PER_CLASS=2, Cd).
- Cu probe runs as a separate entry; SCF failure doesn't abort the batch.
- All outputs present: `selected_structures.json`, init/clean/relaxed CIFs, `dft_vs_symmetry.csv`, ρ–E plot.
- Console prints Spearman ρ_s + converged/discarded counts.
- fjh gated correctly: runs if ffc passes the gate, skips with a warning if not.
- No modifications to existing business files (except allowed compat fixes).
