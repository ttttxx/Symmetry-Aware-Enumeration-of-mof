# DFT Validation Pipeline — Design Spec

**Date:** 2026-06-18 (rev 2)  
**Goal:** Enumerate MOF structures by symmetry class (High/Med/Low), sanitize geometries, run CP2K PBE-D3(BJ) single-point or relaxation, and test whether DFT energy correlates inversely with symmetry index ρ within a fixed composition.

---

## 1. Environment (verified)

| Item | Value |
|------|-------|
| CP2K binary | `/opt/homebrew/bin/cp2k.psmp` (v2026.1, MPI+OMP) |
| CP2K_RUN_CMD | `mpirun -np 4 cp2k.psmp` |
| CP2K_DATA_DIR | `/opt/homebrew/Cellar/cp2k/2026.1_1/share/cp2k/data` |
| GFN1-xTB via CP2K | **Not compiled** (no `xtb` flag in `cp2kflags`) |
| Standalone xTB | `xtb-python` + ASE at `/opt/anaconda3` — PBC confirmed: cell-sensitive energy, periodic FIRE converges |

---

## 2. Composition Strategy (Change 1)

**Problem:** ffc has 7 nodes (4×3c + 3×4c). With two choices per degree, colorings span 20 distinct compositions; mixing them conflates composition effects with symmetry effects.

**Solution:** Enumerate all ffc structures, then restrict to the two compositions that span all three symmetry classes. Analyze each composition independently.

### ffc fixed compositions (both used, never pooled)

| Label | 3c_B_Ch | 3c_triazine_Ch | 4c_Cd_1_Ch | 4c_B_Ch | Total structures | Stab coverage |
|-------|---------|----------------|------------|---------|-----------------|---------------|
| **Comp A** | 2 | 2 | 1 | 2 | 5 | High(3) + Med(1) + Low(1) |
| **Comp B** | 2 | 2 | 2 | 1 | 5 | High(3) + Med(1) + Low(1) |

**Sampling rule for ffc:** Exhaust all 5 structures per composition — no per-class random sampling. `N_PER_CLASS` does **not** apply to ffc.

**Analysis:** Each composition produces its own ρ–E scatter and its own Spearman ρ_s. Results reported as two parallel trends: "Comp A (1 Cd)" and "Comp B (2 Cd)".

### fjh

fjh uses `N_PER_CLASS` per-class random sampling as originally designed (composition is more homogeneous there — to be verified at runtime).

---

## 3. Local vs Server Profiles

### A. Local smoke-test

**Primary run:** `USE_CLOSED_SHELL_VARIANT=True` — swap all `4c_Cu_1_Ch` → `4c_Cd_1_Ch` before enumeration. Cd is closed-shell, uses OT SCF, converges reliably on 16 GB RAM.

**Cu probe:** After the Cd batch, separately select 1 High-sym ffc structure with the **original** Cu `atom_choices` (Cd override bypassed for this entry only). Run as `ENERGY` to exercise the open-shell code path:
- `SPIN_POLARIZED=True`, `SMEARING="auto"` → Fermi-Dirac 300 K, diagonal SCF, UKS
- SCF failure sets `status=failed`; does not abort the batch.
- Config: `CU_PROBE=True`, `CU_PROBE_N=1`

### B. fjh strategy

- `RELAX_MODE = {"ffc": "auto", "fjh": "ENERGY"}` — fjh (~800 atoms) is single-point only.
- `N_PER_CLASS = 2` for fjh locally; bump to 10 on server.
- Gated behind ffc (see §6 analyze.py).

---

## 4. Pre-relaxation Engine

`PRECLEAN_ENGINE = "xtb"` uses standalone **xtb-python** (ASE, PBC-verified):

```python
from xtb.ase.calculator import XTB
from ase.optimize import FIRE

atoms.calc = XTB(method="GFN1-xTB")
opt = FIRE(atoms, logfile=None)
opt.run(fmax=0.1, steps=200)   # 200-step cap prevents divergence from clashing starts
# energy discarded — geometry only
```

Fallback chain: `"xtb"` → `"chgnet"` → `"none"`.

---

## 5. Package Structure

```
dft_validation/
├── __init__.py
├── config.py                # all tunable knobs
├── select_structures.py     # enumerate → filter compositions → assemble CIFs
├── sanitize.py              # wrap → element check → electron count → clash → pre-relax
├── cp2k_runner.py           # generate input → run → parse energy (eV/atom)
├── analyze.py               # post-relax symmetry + per-composition ρ–E + fjh gate
└── run_dft_validation.py    # orchestrator: select→sanitize→run→analyze
```

No existing business files are modified except the two compat fixes in §1 of task.txt, handled in-package.

---

## 6. config.py (final values)

```python
# --- Targets ---
USE_CLOSED_SHELL_VARIANT = True   # True → swap 4c_Cu_1_Ch → 4c_Cd_1_Ch
CU_PROBE      = True              # local only: run 1 Cu High-sym probe for open-shell path test
CU_PROBE_N    = 1

VALIDATION_TARGETS = [
    {"topology": "ffc",
     "atom_choices": {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]},
     "forced_edge": "1B_4H_Ch"},
    {"topology": "fjh",
     "atom_choices": {3: ["3c_B_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]},
     "forced_edge": "1B_4H_Ch"},
]

# ffc: exhaust these two fixed compositions; do NOT apply N_PER_CLASS to ffc
FFC_COMPOSITIONS = [
    {"3c_B_Ch": 2, "3c_triazine_Ch": 2, "4c_Cd_1_Ch": 1, "4c_B_Ch": 2},  # Comp A
    {"3c_B_Ch": 2, "3c_triazine_Ch": 2, "4c_Cd_1_Ch": 2, "4c_B_Ch": 1},  # Comp B
]

# --- Sampling (applies to fjh only) ---
N_PER_CLASS  = 2
SYM_CLASSES  = ["High", "Med", "Low"]
RANDOM_SEED  = 42
MAX_ATOMS    = 250

# --- CHGNet (default off) ---
USE_CHGNET      = False
CHGNET_PRECLEAN = False
CHGNET_COMPARE  = False

# --- CP2K ---
CP2K_RUN_CMD  = "mpirun -np 4 cp2k.psmp"
CP2K_DATA_DIR = "/opt/homebrew/Cellar/cp2k/2026.1_1/share/cp2k/data"

# --- Electronic structure ---
CHARGE         = 0            # assumed electrically neutral framework (confirmed: Cd²⁺ balanced by carboxylates)
SPIN_POLARIZED = "auto"       # auto: metals with unpaired d/f electrons → UKS; else RKS
ODD_ELECTRON_ACTION = "doublet"  # if total electron count is odd: "doublet" (MULTIPLICITY 2) or "flag"
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
FJH_GATE_ON_FFC = True   # qualitative: both comps lowest-E in High AND ρ broadly monotone

# --- Sanitize ---
MIN_DIST_FLOOR    = 0.9
PRECLEAN_ENGINE   = "xtb"   # "xtb" → xtb-python+ASE FIRE | "chgnet" | "none"
SYMPREC_POSTRELAX = 0.3

RESULTS_DIR = "results/dft_validation"
```

---

## 7. Module Responsibilities

### select_structures.py

**ffc path:**
1. Enumerate full structure set using §1 of task.txt pattern (enrich → topo_requirements → MOFEnumerator with degree constraints, `USE_CLOSED_SHELL_VARIANT` applied to atom_choices).
2. For each result, compute its composition key `{ligand: count}`.
3. Retain only structures whose composition matches one of the two entries in `FFC_COMPOSITIONS`. Take all matching structures (no sampling cap).
4. Cu probe: re-enumerate once more with original Cu atom_choices, pick the 1 highest-ρ result from Comp A's 3c pattern, tag as `probe_cu=True`.

**fjh path:**
1. Full enumeration with degree constraints.
2. Bucket by stab value (same logic as `physical_validation.py` lines 105–120).
3. Sample `N_PER_CLASS` per bucket using `RANDOM_SEED`; skip structures with `n_atoms > MAX_ATOMS`.

Both paths: assemble CIFs via `StructureAssembler.assemble_and_save`; save `selected_structures.json` with fields: `topology, composition_label, coloring, stabilizer_size, gamma_size, rho, sym_class, n_atoms, init_cif, probe_cu`.

### sanitize.py

For each init structure, in order:

1. **Wrap + cell check:** pymatgen fractional-coord wrap; verify 3D periodic, non-degenerate cell.
2. **Element coverage check:** every element must have a `DZVP-MOLOPT-SR-GTH` entry in `BASIS_MOLOPT` and a `GTH-PBE` entry in `GTH_POTENTIALS`; missing → `status=missing_kind`, skip.
3. **Electron count parity check:** sum atomic numbers of all atoms minus `CHARGE` = total electrons. If odd → set `multiplicity=2` (doublet), log a warning. If `SPIN_POLARIZED="auto"` hasn't already flagged open-shell, this overrides to UKS + `MULTIPLICITY 2`. Record `n_electrons` and `multiplicity` in the json record.
4. **Minimum-distance check:** `min_dist < MIN_DIST_FLOOR` → pre-relax via `PRECLEAN_ENGINE`. After pre-relax, `min_dist < 0.6 Å` → `status=clash_unfixable`, discard.

Pre-relax via `"xtb"`: `XTB(method="GFN1-xTB")` + `FIRE(fmax=0.1, steps=200)` — periodic, geometry only, energy discarded.

Saves clean CIFs to `<topology>/clean/`; updates json with `status`, `clean_cif`, `n_electrons`, `multiplicity`.

### cp2k_runner.py

Generates CP2K input from clean geometry:

- **Charge/spin:** `CHARGE = config.CHARGE`; multiplicity from sanitize.py record (may be 2 from odd-electron override). `SPIN_POLARIZED` from record: if element list contains Cu/Ni/Fe/Co/Mn/V or multiplicity > 1 → `UKS True`.
- **SCF strategy (auto):**
  - Closed-shell (`UKS False`): OT (`MINIMIZER DIIS`, `PRECONDITIONER FULL_SINGLE_INVERSE`) + `OUTER_SCF`; no smearing.
  - Open-shell/metallic: diagonal diagonalization + `SMEAR METHOD FERMI_DIRAC ELECTRONIC_TEMPERATURE 300` + `ADDED_MOS 20` + `MIXING METHOD BROYDEN_MIXING`; `MAX_SCF 200`, `EPS_SCF 1e-6`.
- `SCF_GUESS ATOMIC` always (no RESTART on first run).
- DFTD3(BJ): `&VDW_POTENTIAL / &PAIR_POTENTIAL TYPE DFTD3 / PARAMETER_FILE_NAME dftd3.dat / REFERENCE_FUNCTIONAL PBE`.
- `RELAX_MODE`: resolved per-topology from config; `"auto"` maps `<100 atoms→CELL_OPT`, `100–300→GEO_OPT`, `>300→ENERGY`. Cu probe always `ENERGY`.
- Γ-point only.
- Runs `CP2K_RUN_CMD` in per-structure subdir; try/except: failure → `status=failed`, log, continue.
- Parses last `ENERGY| Total FORCE_EVAL ... [a.u.]` → energy in Hartree × 27.211386 = eV; `dft_e_per_atom = E_eV / n_atoms`.
- Exports relaxed structure; resumable (skips if `relaxed/*.cif` exists).

### analyze.py

- Re-runs `SpacegroupAnalyzer(symprec=SYMPREC_POSTRELAX)` on relaxed structures.
- **ffc analysis — per composition:**
  - For Comp A and Comp B separately: collect `(rho, dft_e_per_atom)` for converged structures; compute `scipy.stats.spearmanr`; identify lowest-energy structure's sym_class.
  - Produce one ρ–E scatter per composition (sym_class colored), both shown on the same figure as side-by-side subplots.
- **fjh analysis:** single ρ–E scatter, `spearmanr` over all converged.
- **fjh gate (qualitative, evaluated after ffc):** passes if ALL of:
  1. Comp A lowest-energy structure is in the High class.
  2. Comp B lowest-energy structure is in the High class.
  3. Energy is broadly monotone with ρ in both compositions (assessed as: Spearman ρ_s < 0 for both, i.e., higher ρ → lower energy; or at minimum the High class mean < Low class mean for both).
  Returns `(gate_pass: bool, reason: str)` for the orchestrator.
- If `CHGNET_COMPARE`: adds CHGNet-vs-DFT scatter + diagonal + ρ_s annotation.
- Writes `dft_vs_symmetry.csv` (one row per structure): `topology, composition_label, sym_class, rho, dft_e_per_atom, spacegroup_before, spacegroup_after, n_electrons, multiplicity, status`.

### run_dft_validation.py

`main()`:
1. ffc: select → sanitize → run → analyze; prints per-composition Spearman ρ_s + converged/discarded counts.
2. If `FJH_GATE_ON_FFC`: call `analyze.get_fjh_gate()`; print gate result + reason; skip fjh if gate fails.
3. If fjh proceeds: same pipeline with fjh overrides.
4. Terminal summary line per topology: composition label, ρ_s, n_converged/N_total, n_skipped/n_discarded, figure + CSV paths.

---

## 8. Compatibility Fixes (in-package, minimal)

- **NetworkX v3 `shortest_path` returns generator:** monkeypatch before first import of `external_lib` in `select_structures.py`: `import networkx as nx; _sp = nx.shortest_path; nx.shortest_path = lambda *a, **kw: dict(_sp(*a, **kw))`.
- **`enrich_graph_data_with_topology` must precede assembly:** enforced unconditionally in `select_structures.py` (follows `main.py` pattern).

---

## 9. Acceptance

- `python dft_validation/run_dft_validation.py` runs end-to-end on local profile (ffc, 10 Cd structures across Comp A + B, no sampling).
- Cu probe runs as a separate entry; SCF failure doesn't abort the batch.
- All outputs present: `selected_structures.json`, init/clean/relaxed CIFs, `dft_vs_symmetry.csv`, ρ–E figure (2 subplots for ffc).
- Console reports Spearman ρ_s per composition + converged/discarded counts.
- Odd-electron structures flagged in json and run as doublet.
- fjh gated on qualitative two-composition criterion; skips with printed reason if gate fails.
- No modifications to existing business files (except allowed compat fixes).
