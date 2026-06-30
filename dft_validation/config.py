import os

# ════════════════════════════════════════════════════════════════════════════
# SERVER DEPLOYMENT — change ONLY these items on the compute node:
#
#   RUN_DFT_LOCALLY = True          ← enable full GEO_OPT batch
#   START_FROM      = "clean"       ← skip enumeration + sanitize; use synced
#                                      clean/ CIFs + selected_structures.json
#   CP2K_OMP        = <M>           ← OpenMP threads per CP2K job
#   DFT_MAX_WORKERS ← concurrent CP2K jobs per topology
#                                      total active threads ≈ workers * CP2K_OMP
#   MAX_MEM_GB      = <RAM - 4>     ← total memory budget for CP2K jobs
#   CP2K_DATA_DIR   = "/path/to/cp2k/data"   ← server install path
#
# Everything else (cutoffs, basis, XC, RELAX_MODE) stays unchanged.
# On Linux, caffeinate is absent — cp2k_runner detects this automatically.
# ════════════════════════════════════════════════════════════════════════════

# ── Targets ────────────────────────────────────────────────────────────────
USE_CLOSED_SHELL_VARIANT = True   # swap 4c_Cu_1_Ch → 4c_Cd_1_Ch everywhere
CU_PROBE     = True               # one Cu High-sym probe for open-shell path
CU_PROBE_N   = 1

VALIDATION_TARGETS = [
    {
        "topology": "ffc",
        "atom_choices": {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]},
        "forced_edge": "1B_4H_Ch",
    },
    {
        "topology": "fjh",
        "atom_choices": {3: ["3c_B_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]},
        "forced_edge": "1B_4H_Ch",
    },
]

# ffc: exhaust exactly these two compositions (post Cd-substitution names)
FFC_COMPOSITIONS = [
    {"3c_B_Ch": 2, "3c_triazine_Ch": 2, "4c_Cd_1_Ch": 1, "4c_B_Ch": 2},  # Comp A
    {"3c_B_Ch": 2, "3c_triazine_Ch": 2, "4c_Cd_1_Ch": 2, "4c_B_Ch": 1},  # Comp B
]
FFC_COMP_LABELS = ["comp_a", "comp_b"]

# ── Sampling (fjh only — ffc is exhaustive) ────────────────────────────────
N_PER_CLASS  = 2
SYM_CLASSES  = ["High", "Med", "Low"]
RANDOM_SEED  = 42
# ── 结构选择（新口径：组分内 ρ-分层）────────────────────────────────
MAX_ATOMS = 900          # was 250 — 容纳 ~800 原子的 fjh（仅警告，不丢弃）

# ffc：便宜，穷举所有“组分内 ρ≥2”的组分的全部结构
FFC_EXHAUSTIVE = True
FFC_MIN_RHO    = 2       # 只保留组分内有 ≥2 个不同 ρ 的组分
FFC_N_PER_RHO  = 999     # 仅当 FFC_EXHAUSTIVE=False 时生效

# fjh：主验证证据，挑 ρ 最丰富的组分，组分内每个 ρ 抽 N 个
FJH_N_COMPOSITIONS = 2   # 两个固定组分可检验趋势是否可重复
FJH_N_PER_RHO      = 4   # 每个 ρ 值最多抽 4 个；比 2/rho 更稳健
FJH_MIN_RHO        = 3   # 该组分至少要有 3 个不同 ρ 才合格

# 按拓扑 walltime（秒）。fjh 用对角化，O(N³) 比 ffc 慢约两个量级
WALLTIME_SECONDS = {"ffc": 7200, "fjh": 21600}

# ── CHGNet (off by default) ────────────────────────────────────────────────
USE_CHGNET      = False
CHGNET_PRECLEAN = False
CHGNET_COMPARE  = False

# ── CP2K ──────────────────────────────────────────────────────────────────
# RUN_DFT_LOCALLY gates the full 10-structure GEO_OPT batch.
# The local smoke test (1 Cd ENERGY single-point) always runs regardless.
RUN_DFT_LOCALLY = True
# START_FROM: "all" runs the full pipeline (select→sanitize→DFT).
# "clean" skips selection+sanitize and goes straight to DFT using the already-
# synced selected_structures.json and clean/ CIFs — no xtb needed on server.
START_FROM = "clean"
MAX_MEM_GB = 110
CP2K_OMP = 16

# Batch-level concurrency. CP2K itself is already OpenMP-threaded; these
# settings launch multiple independent single-point jobs at once. Keep fjh
# conservative because diagonalization on ~760-atom cells is memory-heavy.
DFT_PARALLEL = True
DFT_MAX_WORKERS = {"ffc": 4, "fjh": 2}
DFT_JOB_MEM_GB = {"ffc": 24, "fjh": 55}
CP2K_DATA_DIR = "/home/huangbo/cp2k/cp2k-2025.2/data"
# SMOKE_CUTOFF_RY is used ONLY for the local single-point smoke test (phase 3).
# Production DFT (phase 4) always uses PW_CUTOFF_RY = 500 Ry.
SMOKE_CUTOFF_RY = 300

# ── Electronic structure ───────────────────────────────────────────────────
CHARGE              = 0            # electrically neutral framework
ODD_ELECTRON_ACTION = "doublet"    # if valence e⁻ count is odd → MULTIPLICITY 2
SPIN_POLARIZED      = "auto"       # auto: Cu/Ni/Fe/Co/Mn/V → UKS; else RKS
DFT_PLUS_U          = None
SMEARING            = "auto"

# ── DFT spec ──────────────────────────────────────────────────────────────
XC_FUNCTIONAL = "PBE"
DISPERSION    = "DFTD3(BJ)"

# Basis: all elements in BASIS_MOLOPT_UZH with pattern DZVP-MOLOPT-PBE-GTH-q{n}
BASIS_SET_FILE = "BASIS_MOLOPT_UZH"
# Per-element: (potential_name, basis_name)
ELEMENT_KIND = {
    "H":  ("GTH-PBE-q1",  "DZVP-MOLOPT-PBE-GTH-q1"),
    "B":  ("GTH-PBE-q3",  "DZVP-MOLOPT-PBE-GTH-q3"),
    "C":  ("GTH-PBE-q4",  "DZVP-MOLOPT-PBE-GTH-q4"),
    "N":  ("GTH-PBE-q5",  "DZVP-MOLOPT-PBE-GTH-q5"),
    "O":  ("GTH-PBE-q6",  "DZVP-MOLOPT-PBE-GTH-q6"),
    "Cu": ("GTH-PBE-q11", "DZVP-MOLOPT-PBE-GTH-q11"),
    "Cd": ("GTH-PBE-q12", "DZVP-MOLOPT-PBE-GTH-q12"),
}
# Valence electrons per element (for parity check)
ELEMENT_VALENCE = {"H": 1, "B": 3, "C": 4, "N": 5, "O": 6, "Cu": 11, "Cd": 12}

PW_CUTOFF_RY = 500
REL_CUTOFF_RY  = 60
# Per-topology relax mode; "auto" → <100 atoms:CELL_OPT, 100-300:GEO_OPT, >300:ENERGY
RELAX_MODE     = {"ffc": "ENERGY", "fjh": "ENERGY"}
HARTREE_TO_EV  = 27.211386

# ── fjh gate ──────────────────────────────────────────────────────────────
FJH_GATE_ON_FFC = True

# ── Sanitize ──────────────────────────────────────────────────────────────
MIN_DIST_FLOOR    = 0.9    # Å — trigger pre-relax if min non-bond dist < this
CLASH_HARD_LIMIT  = 0.6    # Å — discard if still clashing after pre-relax
PRECLEAN_ENGINE   = "xtb"  # "xtb" | "chgnet" | "none"
PRECLEAN_FMAX     = 0.1
PRECLEAN_STEPS    = 200
SYMPREC_POSTRELAX = 0.3

RESULTS_DIR = "results/dft_validation"
