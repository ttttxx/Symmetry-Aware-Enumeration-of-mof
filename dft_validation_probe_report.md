# DFT Validation Probe Report

**Generated:** 2026-06-18  
**Machine:** Darwin 24.6.0 (macOS), 10-core, 16 GB RAM  
**Project root:** `/Users/jing/PyCharmMiscProject/mof`

---

## 1. Environment Summary

### 1.1 CP2K Installation

| Check | Result |
|-------|--------|
| `cp2k` / `cp2k.psmp` / `cp2k.popt` / `cp2k.ssmp` | **NOT FOUND** |
| `module avail cp2k` | `module` command not available |
| `$CP2K_DATA_DIR` | UNSET |

⚠️ **CP2K is not installed on this machine.**

### 1.2 CP2K Data Files

| File | Result |
|------|--------|
| `BASIS_MOLOPT` | NOT FOUND anywhere on disk |
| `GTH_POTENTIALS` | NOT FOUND anywhere on disk |
| `dftd3.dat` | NOT FOUND anywhere on disk |
| `/usr/share/cp2k/data/` | does not exist |
| `/opt/cp2k/data/` | does not exist |
| `/usr/local/cp2k/data/` | does not exist |

### 1.3 Parallel Launcher & Scheduler

| Check | Result |
|-------|--------|
| `mpirun` | `/opt/anaconda3/bin/mpirun` (OpenMPI from conda) |
| `srun` | NOT FOUND |
| `sinfo` / `squeue` / `sbatch` | NOT FOUND |
| SLURM cluster | **None detected** — single-machine only |

### 1.4 Compute Resources

| Metric | Value |
|--------|-------|
| CPU cores (physical) | 10 |
| CPU cores (logical) | 10 |
| Memory | 16 GB (17179869184 bytes) |
| GPU | Not probed (CP2K is CPU-bound) |

---

## 2. Project Structure

### 2.1 Topology Requirements

| Topology | Nodes | Edges | Degree Distribution | \|Aut(G)\| |
|----------|-------|-------|---------------------|------------|
| **ffc** | 7 | 12 | deg=4: 3 nodes, deg=3: 4 nodes | 8 |
| **fjh** | 28 | 48 | deg=4: 12 nodes, deg=3: 16 nodes | 16 |

### 2.2 Available Node CIFs (by degree)

**3c nodes** (for deg-3 positions):

| CIF | Atoms | Key Elements | Metal? |
|-----|-------|-------------|--------|
| `3c_B_Ch.cif` | 9 | C, H | no (organic) |
| `3c_Cu_imadazoline_Ch.cif` | 24 | Cu, N, C, H | **Cu** (3d, open-shell) |
| `3c_pyrrole_Ch.cif` | 7 | N, C, H | no |
| `3c_triazine_Ch.cif` | 6 | C, N | no |
| `3c_uranyl.cif` | 12 | **U**, O, C | **U** (5f, open-shell) |

**4c nodes** (for deg-4 positions):

| CIF | Atoms | Key Elements | Metal? |
|-----|-------|-------------|--------|
| `4c_B_Ch.cif` | 8 | C, H | no (organic) |
| `4c_Cu_1_Ch.cif` | 14 | **Cu**, O, C | **Cu** (3d, open-shell) |
| `4c_Cd_1_Ch.cif` | 13 | **Cd**, O, C | Cd (d10, closed-shell) |
| `4c_Ni_1_Ch.cif` | 14 | **Ni**, O, C | **Ni** (3d, open-shell) |
| `4c_V_1_Ch.cif` | 14 | **V**, O, C | **V** (3d, open-shell) |
| `4c_In_1_Ch.cif` | 13 | **In**, O, C | In (p-block, closed-shell) |
| `4c_Zn_tet.cif` | 5 | Zn, X | Zn (d10, closed-shell) |
| `4c_ethene_Ch.cif` | 2 | C | no |
| `4c_porphyrin_Ch.cif` | 34 | N, C, H | no |
| `4c_pyrene_Ch.cif` | 22 | C, H | no |
| `4c_bicyclooctane_Ch.cif` | 12 | C, H | no |

### 2.3 Available Edge CIFs

| CIF | Atoms | Key Elements |
|-----|-------|-------------|
| `1B_4H_Ch.cif` | 10 | C, H |
| `1B_1TrU.cif` | 12 | C, H |
| `BDC.cif` | 12 | C, H |

### 2.4 Reference from main.py

```python
# main.py line 309-311 (for vba topology)
ATOM_CHOICES = ["4c_B_Ch", "4c_Cd_1_Ch"]
# ATOM_CHOICES = ["4c_Cu_1_Ch", "4c_B_Ch", "3c_B_Ch"]   # commented-out alternative
FORCED_EDGE = "1B_4H_Ch"
```

---

## 3. Proposed Validation Targets

### 3.1 ffc Topology

| Configuration | Value | Rationale |
|--------------|-------|-----------|
| **atom_choices** | `["4c_Cu_1_Ch", "4c_B_Ch", "3c_B_Ch"]` | Cu paddlewheel at deg-4 sites; B-based organic at deg-3; matches commented-out main.py choice |
| **alternate (closed-shell)** | `["4c_Cd_1_Ch", "4c_B_Ch", "3c_B_Ch"]` | Cd replaces Cu — no spin, simpler DFT; main.py's primary choice uses Cd |
| **alternate (all-organic)** | `["4c_B_Ch", "3c_B_Ch"]` | No metals, fastest DFT but less chemically interesting |
| **forced_edge** | `"1B_4H_Ch"` | Consistent with main.py; simple phenyl edge |

⚠️ 化学合理性需人工确认 — see §6.

### 3.2 fjh Topology

| Configuration | Value | Rationale |
|--------------|-------|-----------|
| **atom_choices** | `["4c_Cu_1_Ch", "4c_B_Ch", "3c_B_Ch"]` | Same ligand palette as ffc; 12 × deg-4 + 16 × deg-3 |
| **alternate** | `["4c_Cd_1_Ch", "4c_B_Ch", "3c_B_Ch"]` | Closed-shell variant |
| **forced_edge** | `"1B_4H_Ch"` | Consistent edge choice |

⚠️ fjh is **4× larger** than ffc; DFT may not be feasible on local hardware — see §4.

---

## 4. Assembled Structure Scale

### 4.1 Assembly Results (highest-symmetry coloring)

| Topology | atom_choices | Atoms | Volume (Å³) | Lattice (a,b,c) | Composition | Validation |
|----------|-------------|-------|-------------|-----------------|-------------|------------|
| **ffc** (Cu) | 4c_Cu_1_Ch, 4c_B_Ch, 3c_B_Ch | **198** | 16,564 | 33.32, 33.32, 14.92 | C₁₀₈Cu₆H₆₀O₂₄ | ⚠️ 4 clashes (minor) |
| **ffc** (Cd) | 4c_Cd_1_Ch, 4c_B_Ch, 3c_B_Ch | **195** | 16,564 | 33.32, 33.32, 14.92 | C₁₀₈Cd₃H₆₀O₂₄ | ⚠️ 14 clashes |
| **ffc** (org) | 4c_B_Ch, 3c_B_Ch | **180** | 16,564 | 33.32, 33.32, 14.92 | C₁₁₄H₆₆ | ❌ topology fail |
| **fjh** (Cu) | 4c_Cu_1_Ch, 4c_B_Ch, 3c_B_Ch | **792** | 45,255 | 40.00, 40.00, 28.28 | C₄₃₂Cu₂₄H₂₄₀O₉₆ | ✅ all passed |
| **fjh** (est.) | Cd variant | ~780 | ~45,000 | same | C₄₃₂Cd₁₂H₂₄₀O₉₆ | — not assembled |

*Assembly used `enrich_graph_data_with_topology` (required for `alpha_matrix`/`Bstar`) + NetworkX v3 compatibility fix (`dict(nx.shortest_path(TG))`).*

### 4.2 Scale Classification

| Topology | Atoms | Classification | Feasible on 10-core/16 GB? |
|----------|-------|---------------|---------------------------|
| **ffc** | ~200 | Medium | ⚠️ Borderline — GEO_OPT possible with careful settings |
| **fjh** | ~800 | Large | ❌ Unlikely — would need ENERGY mode or cluster |

---

## 5. DFT Configuration Recommendations

### 5.1 Master Configuration Table

| Parameter | Recommended Value | Reason |
|-----------|------------------|--------|
| **CP2K_RUN_CMD** | `mpirun -np 4 cp2k.psmp` (local) | 10 phys cores; 4 MPI × 2 OMP threads; conservative for 16 GB |
| **CP2K_RUN_CMD** (cluster) | `srun cp2k.psmp` (see §7 sbatch template) | If SLURM cluster becomes available |
| **CP2K_DATA_DIR** | `UNKNOWN` — must be set after install | Required for BASIS_MOLOPT, GTH_POTENTIALS, dftd3.dat |
| **SPIN_POLARIZED** | `True` for Cu/Ni/V/U-containing choices; `False` for Cd/In/Zn/organic-only | Cu²⁺ (d⁹) and Ni²⁺ (d⁸) are open-shell; Cd²⁺ (d¹⁰) is closed-shell |
| **DFT_PLUS_U** | `{"Cu": 4.0}` (example); `None` for closed-shell | Cu 3d electrons are strongly correlated; U=4.0 eV is a common literature value for Cu paddlewheel — ⚠️ needs literature validation |
| **DFT_PLUS_U** (Ni variant) | `{"Ni": 6.0}` (example) | Ni²⁺ often uses U=5–7 eV — ⚠️ needs literature validation |
| **RELAX_MODE** (ffc, ~200 atoms) | **GEO_OPT** (ion relaxation only, fixed cell) | Cell from assembly should be close to equilibrium; ion-only saves ~3–5× vs CELL_OPT |
| **RELAX_MODE** (fjh, ~800 atoms) | **ENERGY** (single-point on CHGNet-relaxed geometry) | Too large for full DFT relaxation on local hardware; CHGNet pre-relaxation validated in physical_validation.py |
| **MAX_ATOMS** (for screening) | 2500 | Same as physical_validation.py default |
| **N_PER_QUADRANT** (ffc) | 3–4 | ~200 atoms, generous sampling affordable |
| **N_PER_QUADRANT** (fjh) | 2 | ~800 atoms, minimize DFT calls |
| **Walltime estimate** (ffc GEO_OPT) | ~4–12 hours per structure (order-of-magnitude) | CP2K GEO_OPT ~200 atoms, 4 MPI ranks, estimated 50–150 SCF steps |
| **Walltime estimate** (fjh ENERGY) | ~1–2 hours per structure | Single-point DFT ~800 atoms |

### 5.2 Relaxation Protocol by Atom Count

```
< 100 atoms  → CELL_OPT (full variable-cell relaxation)
100–300      → GEO_OPT  (ion relaxation, fixed cell)     ← ffc
300–600      → GEO_OPT with tight memory budget, or ENERGY
> 600        → ENERGY  (single-point on pre-relaxed geom)  ← fjh
```

### 5.3 Recommended CP2K Input Fragments

#### For Cu-containing MOF (spin-polarized + DFT+U):

```fortran
&DFT
  BASIS_SET_FILE_NAME  BASIS_MOLOPT
  POTENTIAL_FILE_NAME  GTH_POTENTIALS
  CHARGE 0
  MULTIPLICITY 1              ! will be overridden per metal count
  SPIN_POLARIZED .TRUE.
  &MGRID
    CUTOFF 600
    REL_CUTOFF 60
  &END MGRID
  &QS
    METHOD GPW
    EPS_DEFAULT 1.0E-12
  &END QS
  &SCF
    SCF_GUESS RESTART
    MAX_SCF 200
    EPS_SCF 1.0E-6
    &OT
      MINIMIZER DIIS
      PRECONDITIONER FULL_SINGLE_INVERSE
    &END OT
    &OUTER_SCF
      MAX_SCF 20
      EPS_SCF 1.0E-6
    &END OUTER_SCF
  &END SCF
  &XC
    &XC_FUNCTIONAL PBE
    &END XC_FUNCTIONAL
    &VDW_POTENTIAL
      POTENTIAL_TYPE PAIR_POTENTIAL
      &PAIR_POTENTIAL
        TYPE DFTD3
        PARAMETER_FILE_NAME dftd3.dat
        REFERENCE_FUNCTIONAL PBE
      &END PAIR_POTENTIAL
    &END VDW_POTENTIAL
  &END XC
  &PLUS_U_METHOD
    &ATOMIC_KIND Cu
      L 2            ! d-orbital
      U_MINUS_J 4.0  ! effective U (needs validation)
    &END ATOMIC_KIND
  &END PLUS_U_METHOD
&END DFT
```

---

## 6. Items Requiring Human Review

| # | Item | Why CC Cannot Decide Alone |
|---|------|---------------------------|
| 1 | **atom_choices chemical validity** | The proposed Cu/B combinations for ffc/fjh are based on graph degree compatibility + existing main.py choices. Whether Cu paddlewheel + B-tricarboxylate is chemically synthesizable needs expert judgment. |
| 2 | **DFT+U value for Cu** | U=4.0 eV is a common value for CuO/Cu₂O but the effective U depends on oxidation state and coordination environment. Literature search required (e.g., *Phys. Rev. B* 2015). |
| 3 | **DFT+U for Ni (if used)** | Ni²⁺ U values range 5–8 eV depending on coordination. |
| 4 | **Spin state / MULTIPLICITY** | Cu²⁺ paddlewheel may couple antiferromagnetically. The total multiplicity and spin configuration need manual specification. |
| 5 | **Final MPI rank count** | `mpirun -np 4` is a starting guess for 10-core/16 GB. Memory-per-rank for ~200 atoms with DZVP basis is ~2–3 GB; adjust after first test run. |
| 6 | **fjh feasibility** | 792 atoms on 16 GB with CP2K is very tight. Even ENERGY mode may swap. Consider splitting into smaller sub-systems or waiting for cluster access. |
| 7 | **Assembly quality** | ffc assemblies had minor atom clashes (H–H distances ~1.35 Å). These typically resolve during GEO_OPT but the starting geometry may have high forces. |
| 8 | **Forced edge suitability** | `1B_4H_Ch` (phenyl edge) was validated for vba topology in main.py. Whether it's geometrically suitable for ffc/fjh needs visual inspection of assembled CIFs. |

---

## 7. Missing Environment Components

### CP2K Not Installed — Three Paths Forward

#### Path A: Load Cluster Module (if cluster becomes available)
```bash
module avail cp2k          # find available versions
module load cp2k/2024.1    # example
export CP2K_DATA_DIR=/path/to/cp2k/data
```

#### Path B: Local Installation (macOS)
```bash
# Option 1: Homebrew (simplest)
brew install cp2k

# Option 2: Conda (if brew unavailable)
conda install -c conda-forge cp2k

# Option 3: From source (full control, time-consuming)
# https://github.com/cp2k/cp2k/blob/master/INSTALL.md
# Requires: gfortran, openmpi, fftw, libint, libxc, scalapack, elpa...
```
After install, set:
```bash
export CP2K_DATA_DIR=/usr/local/share/cp2k/data  # or wherever BASIS_MOLOPT lives
```

#### Path C: Apply for HPC Access
- Chinese facilities: 广州超算 (天河), 上海超算, 北京超算
- International: XSEDE (US), PRACE (EU)

### SLURM sbatch Template (if cluster becomes available)

```bash
#!/bin/bash
#SBATCH --job-name=cp2k_mof
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --output=cp2k_%j.out

module load cp2k/2024.1
export CP2K_DATA_DIR=/path/to/cp2k/data

srun cp2k.psmp -i input.inp -o output.out
```

---

## Appendix A: Assembly Code Fixes Applied

Two compatibility fixes were needed to run the assembler with modern Python/NetworkX:

1. **`periodic_structure_analysis.py:874`** — `nx.shortest_path(TG)` → `dict(nx.shortest_path(TG))` (NetworkX v3 returns generator)
2. **`physical_validation.py` must call `enrich_graph_data_with_topology(graph_data)`** to compute `alpha_matrix` and `Bstar` before assembly (otherwise matrix shapes mismatch)

## Appendix B: Probe Script

The assembly probe script is at `_step3_probe.py` (temporary, can be deleted).

## Appendix C: Node/Edge Atom Count Reference

Full per-CIF atom counts used for estimates:

| File | Atoms | Elements |
|------|-------|----------|
| `4c_Cu_1_Ch.cif` | 14 | C:4, Cu:2, O:8 |
| `4c_Cd_1_Ch.cif` | 13 | C:4, Cd:1, O:8 |
| `4c_Ni_1_Ch.cif` | 14 | C:4, Ni:2, O:8 |
| `4c_V_1_Ch.cif` | 14 | C:4, V:2, O:8 |
| `4c_In_1_Ch.cif` | 13 | C:4, In:1, O:8 |
| `4c_B_Ch.cif` | 8 | C:6, H:2 |
| `3c_B_Ch.cif` | 9 | C:6, H:3 |
| `3c_Cu_imadazoline_Ch.cif` | 24 | C:9, Cu:3, H:6, N:6 |
| `3c_uranyl.cif` | 12 | C:3, O:8, U:1 |
| `1B_4H_Ch.cif` | 10 | C:6, H:4 |
