import os

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
MAX_ATOMS    = 250

# ── CHGNet (off by default) ────────────────────────────────────────────────
USE_CHGNET      = False
CHGNET_PRECLEAN = False
CHGNET_COMPARE  = False

# ── CP2K ──────────────────────────────────────────────────────────────────
CP2K_RUN_CMD  = "mpirun -np 4 cp2k.psmp"
CP2K_DATA_DIR = "/opt/homebrew/Cellar/cp2k/2026.1_1/share/cp2k/data"

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

PW_CUTOFF_RY   = 600
REL_CUTOFF_RY  = 60
# Per-topology relax mode; "auto" → <100 atoms:CELL_OPT, 100-300:GEO_OPT, >300:ENERGY
RELAX_MODE     = {"ffc": "auto", "fjh": "ENERGY"}
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
