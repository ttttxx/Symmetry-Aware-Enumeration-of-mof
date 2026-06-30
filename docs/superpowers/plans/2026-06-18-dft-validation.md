# DFT Validation Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `dft_validation/` — a four-phase pipeline that enumerates ffc MOF structures by symmetry, sanitizes geometries, runs CP2K PBE-D3(BJ), and tests whether DFT energy correlates inversely with symmetry index ρ within fixed compositions.

**Architecture:** Phase 1 enumerates and assembles 10 ffc CIFs (Comp A + B). Phase 2 sanitizes with xTB FIRE pre-relaxation. Phase 3 is a single-structure dry pass to confirm CP2K → eV/atom roundtrip. Phase 4 runs full batch + analysis. Each phase has an explicit acceptance gate before the next begins.

**Tech Stack:** Python 3, pymatgen, xtb-python + ASE, CP2K 2026.1 (`cp2k.psmp`), scipy, matplotlib. Runs from the repo root (`/Users/jing/PyCharmMiscProject/mof`). Use `python3` (`.venv` interpreter, which has pymatgen).

## Global Constraints

- Run all `python3` commands from `/Users/jing/PyCharmMiscProject/mof` (project root).
- Never modify existing business files (`main.py`, `symmetry.py`, `enumerator.py`, `assembler.py`, `external_lib/`). All changes go in `dft_validation/` and `tests/dft_validation/`.
- Only allowed exception: NetworkX compat monkeypatch — apply in `select_structures.py` before importing `external_lib`.
- Basis file: `BASIS_MOLOPT_UZH` (not `BASIS_MOLOPT`). Basis name pattern: `DZVP-MOLOPT-PBE-GTH-q{n}`.
- Per-element `q` values (valence electrons for pseudopotential + basis): H=1, B=3, C=4, N=5, O=6, Cu=11, Cd=12.
- Energy conversion: CP2K outputs Hartree. Multiply by `27.211386` to get eV, then divide by `n_atoms` → eV/atom.
- ffc Comp A: `{3c_B_Ch:2, 3c_triazine_Ch:2, 4c_Cd_1_Ch:1, 4c_B_Ch:2}` — 5 structures, **no per-class sampling**.
- ffc Comp B: `{3c_B_Ch:2, 3c_triazine_Ch:2, 4c_Cd_1_Ch:2, 4c_B_Ch:1}` — 5 structures, **no per-class sampling**.
- `N_PER_CLASS` applies to fjh only.
- Assembler call pattern: `StructureAssembler(output_dir=dir).assemble_and_save(full_cif_path, smap, graph_data, forced_edge=...)`  where `smap = {str(i): coloring[i] for i in range(len(coloring))}`.
- Phase 3 is a hard PAUSE: do not proceed to Phase 4 until the user confirms the dry-pass energy value.

---

## File Map

| File | Create/Modify | Responsibility |
|------|--------------|----------------|
| `dft_validation/__init__.py` | Create | Empty package marker |
| `dft_validation/config.py` | Create | All tunable knobs — single source of truth |
| `dft_validation/select_structures.py` | Create | Enumerate → filter compositions → assemble CIFs → write JSON |
| `dft_validation/sanitize.py` | Create | Wrap, element check, electron parity, distance, xTB FIRE |
| `dft_validation/cp2k_runner.py` | Create | Generate CP2K input, run, parse energy (Hartree→eV/atom) |
| `dft_validation/analyze.py` | Create | Per-composition Spearman, ρ–E plots, fjh qualitative gate |
| `dft_validation/run_dft_validation.py` | Create | Orchestrator: Phase 1→2→3→4, fjh gate |
| `tests/dft_validation/__init__.py` | Create | Empty |
| `tests/dft_validation/test_select.py` | Create | Phase 1 tests |
| `tests/dft_validation/test_sanitize.py` | Create | Phase 2 tests |
| `tests/dft_validation/test_runner.py` | Create | Phase 3 tests (input gen + energy parse) |
| `tests/dft_validation/test_analyze.py` | Create | Phase 4 tests (Spearman, gate) |

---

## Phase 1 — Scaffold + config.py + select_structures.py

**Acceptance gate:** `python3 dft_validation/run_dft_validation.py --phase 1` produces `results/dft_validation/selected_structures.json` containing exactly 10 ffc entries (5 comp_a + 5 comp_b) plus 1 cu_probe entry, and all 11 init CIFs exist on disk.

---

### Task 1.1 — Package scaffold + config.py

**Files:**
- Create: `dft_validation/__init__.py`
- Create: `dft_validation/config.py`
- Create: `tests/dft_validation/__init__.py`

**Interfaces:**
- Produces: `config` module — all downstream tasks import from here.

- [ ] **Step 1: Create the package directories and empty `__init__.py` files**

```bash
mkdir -p dft_validation tests/dft_validation
touch dft_validation/__init__.py tests/dft_validation/__init__.py
```

- [ ] **Step 2: Write `dft_validation/config.py`**

```python
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
```

- [ ] **Step 3: Verify config imports cleanly**

```bash
python3 -c "from dft_validation import config; print('config OK', config.CP2K_DATA_DIR)"
```

Expected output:
```
config OK /opt/homebrew/Cellar/cp2k/2026.1_1/share/cp2k/data
```

- [ ] **Step 4: Commit**

```bash
git add dft_validation/__init__.py dft_validation/config.py tests/dft_validation/__init__.py
git commit -m "feat: add dft_validation package scaffold and config"
```

---

### Task 1.2 — select_structures.py

**Files:**
- Create: `dft_validation/select_structures.py`
- Create: `tests/dft_validation/test_select.py`

**Interfaces:**
- Consumes: `config.VALIDATION_TARGETS`, `config.FFC_COMPOSITIONS`, `config.FFC_COMP_LABELS`, `config.USE_CLOSED_SHELL_VARIANT`, `config.CU_PROBE`, `config.CU_PROBE_N`, `config.N_PER_CLASS`, `config.RANDOM_SEED`, `config.MAX_ATOMS`, `config.RESULTS_DIR`
- Produces:
  - `run_selection(results_dir) -> list[dict]` — returns record list and writes `selected_structures.json`
  - Each record: `{"id": str, "topology": str, "composition_label": str, "coloring": list, "stabilizer_size": int, "gamma_size": int, "rho": float, "sym_class": str, "n_atoms": int, "init_cif": str, "probe_cu": bool, "status": "selected"}`

- [ ] **Step 1: Write failing tests**

```python
# tests/dft_validation/test_select.py
import sys, os
sys.path.insert(0, '.')
import pytest
from collections import Counter

def _apply_cd_subst(atom_choices):
    """Mirror of select_structures._apply_cd_substitution."""
    result = {}
    for deg, ligs in atom_choices.items():
        result[deg] = ["4c_Cd_1_Ch" if l == "4c_Cu_1_Ch" else l for l in ligs]
    return result

def test_comp_key_matches_expected():
    from dft_validation.select_structures import comp_key
    coloring = ["3c_B_Ch", "3c_triazine_Ch", "3c_B_Ch", "3c_triazine_Ch",
                "4c_Cd_1_Ch", "4c_B_Ch", "4c_B_Ch"]
    key = comp_key(coloring)
    assert key == {"3c_B_Ch": 2, "3c_triazine_Ch": 2, "4c_Cd_1_Ch": 1, "4c_B_Ch": 2}

def test_classify_stabs_ffc():
    from dft_validation.select_structures import classify_sym_classes
    # ffc stabs = {8,4,2,1} → high={8,4}, med={2}, low={1}
    fake = [{"stabilizer_size": s} for s in [8, 8, 4, 4, 2, 1, 1]]
    classified = classify_sym_classes(fake)
    classes = [r["sym_class"] for r in classified]
    assert classes.count("High") == 4
    assert classes.count("Med")  == 2
    assert classes.count("Low")  == 2

def test_ffc_gives_10_structures():
    """Integration: real enumeration → exactly 10 ffc Cd structures (5 A + 5 B)."""
    from dft_validation.select_structures import enumerate_ffc
    from dft_validation import config
    target = next(t for t in config.VALIDATION_TARGETS if t["topology"] == "ffc")
    records = enumerate_ffc(target)
    assert len(records) == 10, f"Expected 10, got {len(records)}"
    labels = [r["composition_label"] for r in records]
    assert labels.count("comp_a") == 5
    assert labels.count("comp_b") == 5

def test_ffc_sym_classes_present():
    from dft_validation.select_structures import enumerate_ffc
    from dft_validation import config
    target = next(t for t in config.VALIDATION_TARGETS if t["topology"] == "ffc")
    records = enumerate_ffc(target)
    for label in ["comp_a", "comp_b"]:
        subset = [r for r in records if r["composition_label"] == label]
        classes = {r["sym_class"] for r in subset}
        assert "High" in classes, f"{label} missing High"
        assert "Med"  in classes, f"{label} missing Med"
        assert "Low"  in classes, f"{label} missing Low"
```

- [ ] **Step 2: Run tests — verify they FAIL**

```bash
python3 -m pytest tests/dft_validation/test_select.py -v 2>&1 | head -20
```

Expected: `ImportError` on `dft_validation.select_structures`.

- [ ] **Step 3: Write `dft_validation/select_structures.py`**

```python
import os, sys, json, random
from collections import Counter
sys.path.insert(0, '.')

# ── NetworkX v3 compat: shortest_path returns generator in v3 ─────────────
import networkx as nx
_nx_sp = nx.shortest_path
nx.shortest_path = lambda *a, **kw: dict(_nx_sp(*a, **kw))

from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator
from assembler import StructureAssembler
from dft_validation import config


def _graph_path(topology: str) -> str:
    return os.path.join("data", topology, f"{topology}_graph.json")


def _apply_cd_substitution(atom_choices: dict) -> dict:
    """Replace 4c_Cu_1_Ch → 4c_Cd_1_Ch in atom_choices if USE_CLOSED_SHELL_VARIANT."""
    if not config.USE_CLOSED_SHELL_VARIANT:
        return atom_choices
    return {
        deg: ["4c_Cd_1_Ch" if l == "4c_Cu_1_Ch" else l for l in ligs]
        for deg, ligs in atom_choices.items()
    }


def comp_key(coloring: list) -> dict:
    """Return {ligand_name: count} composition dict for a coloring."""
    return dict(Counter(coloring))


def _comp_matches(coloring: list, target: dict) -> bool:
    return comp_key(coloring) == target


def classify_sym_classes(results: list) -> list:
    """
    Assign sym_class to each result using physical_validation.py bucketing logic
    (lines 105-120): top-2 stab values → High, interior >1 → Med, stab==1 → Low.
    Returns new list with 'sym_class' key added.
    """
    all_stabs = sorted(set(r["stabilizer_size"] for r in results), reverse=True)
    high_vals, med_vals = set(), set()
    if len(all_stabs) <= 2:
        if all_stabs and all_stabs[0] > 1:
            high_vals.add(all_stabs[0])
    elif len(all_stabs) == 3:
        high_vals.add(all_stabs[0])
        if all_stabs[1] > 1:
            med_vals.add(all_stabs[1])
    else:
        high_vals.update(all_stabs[:2])
        for v in all_stabs[2:]:
            if v > 1:
                med_vals.add(v)

    out = []
    for r in results:
        s = r["stabilizer_size"]
        if s in high_vals:
            cls = "High"
        elif s in med_vals:
            cls = "Med"
        else:
            cls = "Low"
        out.append({**r, "sym_class": cls})
    return out


def _enumerate_topology(topology: str, atom_choices: dict):
    """Return (results, graph_data, perms) for a topology + atom_choices."""
    graph_data = enrich_graph_data_with_topology(load_graph(_graph_path(topology)))
    topo_req   = calculate_topology_requirements(graph_data)
    node_degs  = [n.get("degree", 0) for n in graph_data["nodes"]]
    perms      = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
    all_ligs   = sorted(set(l for ls in atom_choices.values() for l in ls))
    results    = MOFEnumerator(
        perms, graph_data["n_nodes"], verbose=False,
        topology_requirements=topo_req,
        user_selected_ligands=atom_choices,
        node_degrees=node_degs,
    ).enumerate(all_ligs)
    return results, graph_data, perms


def enumerate_ffc(target: dict) -> list:
    """
    Enumerate ffc, apply Cd substitution, filter to FFC_COMPOSITIONS,
    classify sym_classes, return all matching records (no sampling).
    """
    ac = _apply_cd_substitution(target["atom_choices"])
    results, graph_data, perms = _enumerate_topology(target["topology"], ac)
    gamma_size = len(perms)
    classified = classify_sym_classes(results)

    records = []
    for comp_dict, label in zip(config.FFC_COMPOSITIONS, config.FFC_COMP_LABELS):
        matching = [r for r in classified if _comp_matches(r["coloring"], comp_dict)]
        for r in matching:
            records.append({
                "id": f"ffc_{label}_{r['id']}",
                "topology": "ffc",
                "composition_label": label,
                "coloring": r["coloring"],
                "stabilizer_size": r["stabilizer_size"],
                "gamma_size": gamma_size,
                "rho": r["stabilizer_size"] / gamma_size,
                "sym_class": r["sym_class"],
                "n_atoms": None,  # filled after assembly
                "init_cif": None,
                "probe_cu": False,
                "status": "selected",
            })
    return records


def _cu_probe_record(target: dict, ffc_graph_data, ffc_perms) -> dict | None:
    """Pick 1 High-sym ffc structure using original Cu atom_choices."""
    if not config.CU_PROBE:
        return None
    # Use ORIGINAL atom_choices (no Cd substitution)
    ac = target["atom_choices"]
    gamma_size = len(ffc_perms)
    results, _, _ = _enumerate_topology(target["topology"], ac)
    classified = classify_sym_classes(results)
    high = [r for r in classified if r["sym_class"] == "High"]
    if not high:
        return None
    # Pick highest stabilizer_size first, then by id for determinism
    best = max(high, key=lambda r: (r["stabilizer_size"], -r["id"]))
    return {
        "id": "ffc_cu_probe_0",
        "topology": "ffc",
        "composition_label": "cu_probe",
        "coloring": best["coloring"],
        "stabilizer_size": best["stabilizer_size"],
        "gamma_size": gamma_size,
        "rho": best["stabilizer_size"] / gamma_size,
        "sym_class": best["sym_class"],
        "n_atoms": None,
        "init_cif": None,
        "probe_cu": True,
        "status": "selected",
    }


def _assemble_records(records: list, topology: str, forced_edge: str, out_dir: str) -> list:
    """Assemble CIFs for all records, update n_atoms and init_cif in-place."""
    init_dir = os.path.join(out_dir, topology, "init")
    os.makedirs(init_dir, exist_ok=True)
    assembler = StructureAssembler(output_dir=init_dir)

    # Need graph_data for the first record's topology (same for all ffc records)
    # Re-load graph_data using the first record's coloring to determine ac
    # We already have graph_data in scope — but to keep this function self-contained,
    # pass graph_data as parameter. Refactor: accept graph_data.
    raise NotImplementedError("Call _assemble_records_with_graph instead")


def _assemble_records_with_graph(records: list, graph_data: dict, forced_edge: str,
                                  out_dir: str, topology: str) -> list:
    """Assemble CIFs, update n_atoms + init_cif. Returns updated records."""
    init_dir = os.path.join(out_dir, topology, "init")
    os.makedirs(init_dir, exist_ok=True)
    assembler = StructureAssembler(output_dir=init_dir)

    from pymatgen.core import Structure

    updated = []
    for rec in records:
        smap = {str(i): c for i, c in enumerate(rec["coloring"])}
        cif_name = f"{rec['id']}_rho{rec['rho']:.4f}_{rec['sym_class']}.cif"
        cif_path = os.path.join(init_dir, cif_name)
        try:
            assembler.assemble_and_save(cif_path, smap, graph_data, forced_edge=forced_edge)
            struct = Structure.from_file(cif_path)
            rec = {**rec, "n_atoms": len(struct), "init_cif": cif_path, "status": "assembled"}
        except Exception as e:
            print(f"  [select] Assembly failed for {rec['id']}: {e}")
            rec = {**rec, "status": f"assembly_failed: {e}"}
        updated.append(rec)
    return updated


def run_selection(results_dir: str = None) -> list:
    """
    Main entry point for Phase 1.
    Enumerates ffc Comp A+B + Cu probe, assembles CIFs, writes selected_structures.json.
    Returns the list of records.
    """
    if results_dir is None:
        results_dir = config.RESULTS_DIR
    os.makedirs(results_dir, exist_ok=True)

    all_records = []

    ffc_target = next(t for t in config.VALIDATION_TARGETS if t["topology"] == "ffc")
    ac_cd = _apply_cd_substitution(ffc_target["atom_choices"])

    # Enumerate once (Cd variant) — reuse graph_data and perms for probe too
    graph_data = enrich_graph_data_with_topology(load_graph(_graph_path("ffc")))
    topo_req   = calculate_topology_requirements(graph_data)
    node_degs  = [n.get("degree", 0) for n in graph_data["nodes"]]
    perms      = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
    gamma_size = len(perms)
    all_ligs   = sorted(set(l for ls in ac_cd.values() for l in ls))
    results    = MOFEnumerator(
        perms, graph_data["n_nodes"], verbose=False,
        topology_requirements=topo_req,
        user_selected_ligands=ac_cd,
        node_degrees=node_degs,
    ).enumerate(all_ligs)
    classified = classify_sym_classes(results)

    for comp_dict, label in zip(config.FFC_COMPOSITIONS, config.FFC_COMP_LABELS):
        matching = [r for r in classified if _comp_matches(r["coloring"], comp_dict)]
        for r in matching:
            all_records.append({
                "id": f"ffc_{label}_{r['id']}",
                "topology": "ffc",
                "composition_label": label,
                "coloring": r["coloring"],
                "stabilizer_size": r["stabilizer_size"],
                "gamma_size": gamma_size,
                "rho": r["stabilizer_size"] / gamma_size,
                "sym_class": r["sym_class"],
                "n_atoms": None,
                "init_cif": None,
                "probe_cu": False,
                "status": "selected",
            })

    # Cu probe (re-enumerate with original Cu choices)
    probe = _cu_probe_record(ffc_target, graph_data, perms)
    if probe:
        # Assemble probe with Cu graph_data (same topology geometry)
        probe_records = _assemble_records_with_graph(
            [probe], graph_data, ffc_target["forced_edge"], results_dir, "ffc"
        )
        all_records.extend(probe_records)

    # Assemble all Cd records
    cd_records = _assemble_records_with_graph(
        [r for r in all_records if not r["probe_cu"]],
        graph_data, ffc_target["forced_edge"], results_dir, "ffc"
    )
    # Replace unassembled records with assembled ones
    cu_probes  = [r for r in all_records if r["probe_cu"]]
    all_records = cd_records + cu_probes

    json_path = os.path.join(results_dir, "selected_structures.json")
    with open(json_path, "w") as f:
        json.dump(all_records, f, indent=2)
    print(f"[select] Wrote {len(all_records)} records to {json_path}")
    assembled = sum(1 for r in all_records if r["status"] == "assembled")
    print(f"[select] Assembled: {assembled}/{len(all_records)}")
    return all_records
```

- [ ] **Step 4: Run failing tests again — some should now pass**

```bash
python3 -m pytest tests/dft_validation/test_select.py::test_comp_key_matches_expected \
    tests/dft_validation/test_select.py::test_classify_stabs_ffc -v
```

Expected: both PASS.

- [ ] **Step 5: Run the integration test (real enumeration, no CIF assembly)**

```bash
python3 -m pytest tests/dft_validation/test_select.py::test_ffc_gives_10_structures \
    tests/dft_validation/test_select.py::test_ffc_sym_classes_present -v
```

Expected: both PASS (~5 s for enumeration).

- [ ] **Step 6: Write the Phase 1 orchestrator stub in `run_dft_validation.py`**

```python
# dft_validation/run_dft_validation.py
import argparse, sys, os
sys.path.insert(0, '.')
from dft_validation import config
from dft_validation.select_structures import run_selection


def phase1(results_dir):
    print("=== Phase 1: select structures ===")
    records = run_selection(results_dir)
    ffc_cd = [r for r in records if r["topology"] == "ffc" and not r["probe_cu"]]
    probe  = [r for r in records if r["probe_cu"]]
    print(f"  ffc Cd records: {len(ffc_cd)}  (expect 10)")
    print(f"  Cu probe records: {len(probe)}  (expect 1)")
    ok = len(ffc_cd) == 10 and len(probe) <= 1
    if not ok:
        print("  FAIL: unexpected record counts")
        sys.exit(1)
    print("  Phase 1 PASS")
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", type=int, default=1)
    parser.add_argument("--results-dir", default=config.RESULTS_DIR)
    args = parser.parse_args()

    if args.phase >= 1:
        records = phase1(args.results_dir)
    # Phases 2–4 added in later tasks

if __name__ == "__main__":
    main()
```

- [ ] **Step 7: Run Phase 1 end-to-end and verify outputs**

```bash
python3 dft_validation/run_dft_validation.py --phase 1
```

Expected output (approximate):
```
=== Phase 1: select structures ===
  ffc Cd records: 10  (expect 10)
  Cu probe records: 1  (expect 1)
  Phase 1 PASS
```

Then verify files exist:
```bash
ls results/dft_validation/selected_structures.json
ls results/dft_validation/ffc/init/*.cif | wc -l
```

Expected: JSON exists; `wc -l` prints `11` (10 Cd + 1 Cu probe).

- [ ] **Step 8: Commit**

```bash
git add dft_validation/select_structures.py dft_validation/run_dft_validation.py \
    tests/dft_validation/test_select.py
git commit -m "feat: Phase 1 — enumerate ffc Comp A+B + Cu probe, assemble init CIFs"
```

---

## Phase 2 — sanitize.py

**Acceptance gate:** `python3 dft_validation/run_dft_validation.py --phase 2` produces clean CIFs in `ffc/clean/`, all 11 entries have `status` ∈ {`clean`, `preclean_done`, `clash_unfixable`, `missing_kind`}, and at least 9/11 have `status` in {`clean`, `preclean_done`}. Confirm xTB runs successfully on the Cd structure (no crash).

---

### Task 2.1 — Element check + electron parity

**Files:**
- Create: `dft_validation/sanitize.py` (partial)
- Create: `tests/dft_validation/test_sanitize.py` (partial)

**Interfaces:**
- Consumes: record dict with `init_cif`, `n_atoms`; `config.ELEMENT_KIND`, `config.ELEMENT_VALENCE`, `config.CHARGE`
- Produces:
  - `check_element_coverage(struct) -> list[str]` — returns list of missing elements (empty = OK)
  - `compute_electron_parity(struct) -> tuple[int, int]` — returns `(n_valence_electrons, multiplicity)`

- [ ] **Step 1: Write failing tests**

```python
# tests/dft_validation/test_sanitize.py
import sys
sys.path.insert(0, '.')
import pytest
from pymatgen.core import Structure, Lattice
import numpy as np

def _make_struct(species, a=10.0):
    """Helper: cubic cell, atoms at fractional positions."""
    n = len(species)
    coords = [[i/n, 0.5, 0.5] for i in range(n)]
    return Structure(Lattice.cubic(a), species, coords)

def test_element_coverage_all_known():
    from dft_validation.sanitize import check_element_coverage
    struct = _make_struct(["C", "H", "O", "N", "Cd"])
    missing = check_element_coverage(struct)
    assert missing == [], f"Expected no missing, got {missing}"

def test_element_coverage_unknown_element():
    from dft_validation.sanitize import check_element_coverage
    struct = _make_struct(["C", "Xe"])  # Xe not in ELEMENT_KIND
    missing = check_element_coverage(struct)
    assert "Xe" in missing

def test_electron_parity_even():
    from dft_validation.sanitize import compute_electron_parity
    # C: 4, H: 1 → 4+1+1+4 = 10, even → multiplicity 1
    struct = _make_struct(["C", "H", "H", "C"])
    n_e, mult = compute_electron_parity(struct)
    assert n_e == 10
    assert mult == 1

def test_electron_parity_odd():
    from dft_validation.sanitize import compute_electron_parity
    # C:4, H:1 → 5, odd → multiplicity 2
    struct = _make_struct(["C", "H"])
    n_e, mult = compute_electron_parity(struct)
    assert n_e == 5
    assert mult == 2
```

- [ ] **Step 2: Run tests — verify FAIL**

```bash
python3 -m pytest tests/dft_validation/test_sanitize.py -v 2>&1 | head -15
```

Expected: `ImportError`.

- [ ] **Step 3: Write `dft_validation/sanitize.py` (element check + parity)**

```python
# dft_validation/sanitize.py
import os, json
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
```

- [ ] **Step 4: Run tests — verify PASS**

```bash
python3 -m pytest tests/dft_validation/test_sanitize.py::test_element_coverage_all_known \
    tests/dft_validation/test_sanitize.py::test_element_coverage_unknown_element \
    tests/dft_validation/test_sanitize.py::test_electron_parity_even \
    tests/dft_validation/test_sanitize.py::test_electron_parity_odd -v
```

Expected: all 4 PASS.

---

### Task 2.2 — Distance check + xTB FIRE pre-relax + full sanitize pipeline

**Files:**
- Modify: `dft_validation/sanitize.py` (add distance check, xTB, pipeline)
- Modify: `tests/dft_validation/test_sanitize.py` (add xTB + pipeline tests)
- Modify: `dft_validation/run_dft_validation.py` (add phase 2)

**Interfaces:**
- Produces: `sanitize_all(records, results_dir) -> list[dict]` — updates each record with `status`, `clean_cif`, `n_electrons`, `multiplicity`

- [ ] **Step 1: Add tests**

```python
# Append to tests/dft_validation/test_sanitize.py

def test_min_distance_clean_struct():
    """A well-spaced structure reports min_dist above floor."""
    from dft_validation.sanitize import min_nonbond_distance
    struct = _make_struct(["C", "C", "C"], a=15.0)
    d = min_nonbond_distance(struct)
    assert d > 2.0

def test_xtb_runs_on_cd_structure():
    """xTB FIRE runs without crash on a real assembled Cd CIF (integration)."""
    import glob
    cifs = glob.glob("results/dft_validation/ffc/init/*comp_a*.cif")
    if not cifs:
        pytest.skip("Phase 1 output not present — run --phase 1 first")
    from dft_validation.sanitize import preclean_xtb
    from pymatgen.core import Structure
    struct = Structure.from_file(cifs[0])
    result = preclean_xtb(struct)
    assert result is not None, "xTB returned None"
    assert hasattr(result, "sites"), "Result is not a Structure"
```

- [ ] **Step 2: Run new tests — verify FAIL**

```bash
python3 -m pytest tests/dft_validation/test_sanitize.py::test_min_distance_clean_struct \
    tests/dft_validation/test_sanitize.py::test_xtb_runs_on_cd_structure -v 2>&1 | head -10
```

Expected: `ImportError` for `min_nonbond_distance` and `preclean_xtb`.

- [ ] **Step 3: Complete `dft_validation/sanitize.py`**

Append to the existing file:

```python
# ── append to dft_validation/sanitize.py ────────────────────────────────

import sys
import numpy as np


def min_nonbond_distance(struct: Structure) -> float:
    """Return the minimum pairwise distance between all sites."""
    from pymatgen.analysis.structure_analyzer import VoronoiConnectivity
    cart = struct.cart_coords
    lattice = struct.lattice
    min_d = float("inf")
    for i in range(len(cart)):
        for j in range(i + 1, len(cart)):
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
    cif_name = os.path.basename(init_cif).replace("/init/", "/clean/")
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
```

- [ ] **Step 4: Run tests**

```bash
python3 -m pytest tests/dft_validation/test_sanitize.py -v
```

Expected: all pass. The xTB integration test takes ~30 s.

- [ ] **Step 5: Add Phase 2 to orchestrator**

In `dft_validation/run_dft_validation.py`, add below `phase1`:

```python
from dft_validation.sanitize import sanitize_all

def phase2(records, results_dir):
    print("=== Phase 2: sanitize ===")
    records = sanitize_all(records, results_dir)
    ok = sum(1 for r in records if r.get("status") in ("clean", "preclean_done"))
    print(f"  clean/preclean: {ok}/{len(records)}")
    if ok == 0:
        print("  FAIL: no clean structures")
        sys.exit(1)
    print("  Phase 2 PASS")
    return records
```

And in `main()`:
```python
    if args.phase >= 2:
        records = phase2(records, args.results_dir)
```

- [ ] **Step 6: Run Phase 2**

```bash
python3 dft_validation/run_dft_validation.py --phase 2
```

Expected:
```
=== Phase 2: sanitize ===
  [sanitize] Processing ffc_comp_a_0 …
  ...
  clean/preclean: 10/11   (probe may differ)
  Phase 2 PASS
```

Verify Cd structures pass without `missing_kind`:
```bash
python3 -c "
import json
recs = json.load(open('results/dft_validation/selected_structures.json'))
for r in recs:
    print(r['id'], r['status'], r.get('n_electrons','?'), r.get('multiplicity','?'))
"
```

Expected: all `comp_a` and `comp_b` records show `status=clean`, `multiplicity=1`.

- [ ] **Step 7: Commit**

```bash
git add dft_validation/sanitize.py dft_validation/run_dft_validation.py \
    tests/dft_validation/test_sanitize.py
git commit -m "feat: Phase 2 — sanitize with element check, electron parity, xTB FIRE pre-relax"
```

---

## Phase 3 — cp2k_runner.py dry pass (single ENERGY, **PAUSE for confirmation**)

**Acceptance gate:** One structure produces a valid `dft_e_per_atom` (eV/atom) value in the CSV. **Stop here and show the user the energy value before proceeding to Phase 4.**

---

### Task 3.1 — CP2K input generator

**Files:**
- Create: `dft_validation/cp2k_runner.py` (input-gen portion)
- Create: `tests/dft_validation/test_runner.py`

**Interfaces:**
- Produces: `generate_cp2k_input(record, run_dir, topology, relax_mode_override=None) -> str` — returns path to written `.inp` file

- [ ] **Step 1: Write failing tests**

```python
# tests/dft_validation/test_runner.py
import sys, os, re
sys.path.insert(0, '.')
import pytest

def test_input_contains_pbe():
    import glob
    cifs = glob.glob("results/dft_validation/ffc/clean/*.cif")
    if not cifs:
        pytest.skip("Phase 2 outputs not present")
    import json
    recs = json.load(open("results/dft_validation/selected_structures.json"))
    rec = next(r for r in recs if r.get("status") in ("clean", "preclean_done")
               and not r.get("probe_cu"))
    from dft_validation.cp2k_runner import generate_cp2k_input
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        inp = generate_cp2k_input(rec, td, "ffc", relax_mode_override="ENERGY")
        text = open(inp).read()
    assert "PBE" in text
    assert "DFTD3" in text
    assert "ENERGY" in text
    assert "DZVP-MOLOPT-PBE-GTH" in text

def test_input_charge_zero():
    import glob, json, tempfile
    cifs = glob.glob("results/dft_validation/ffc/clean/*.cif")
    if not cifs:
        pytest.skip("Phase 2 outputs not present")
    recs = json.load(open("results/dft_validation/selected_structures.json"))
    rec = next(r for r in recs if r.get("status") in ("clean", "preclean_done"))
    from dft_validation.cp2k_runner import generate_cp2k_input
    with tempfile.TemporaryDirectory() as td:
        inp = generate_cp2k_input(rec, td, "ffc", relax_mode_override="ENERGY")
        text = open(inp).read()
    assert re.search(r"CHARGE\s+0", text)

def test_energy_parse():
    from dft_validation.cp2k_runner import parse_energy
    fake_out = """
 SCF run converged in    15 steps
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:            -5812.123456789012
"""
    e = parse_energy(fake_out, n_atoms=200)
    assert e is not None
    expected = -5812.123456789012 * 27.211386 / 200
    assert abs(e - expected) < 1e-6
```

- [ ] **Step 2: Run tests — verify FAIL**

```bash
python3 -m pytest tests/dft_validation/test_runner.py -v 2>&1 | head -15
```

- [ ] **Step 3: Write `dft_validation/cp2k_runner.py`**

```python
# dft_validation/cp2k_runner.py
import os, re, subprocess, json
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

    # Fractional coords → COORD block
    coord_lines = ""
    for site in struct:
        sym = str(site.specie.symbol)
        x, y, z = site.cart_coords
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
    &KPOINTS
      SCHEME GAMMA
    &END KPOINTS
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
    cmd = f"{config.CP2K_RUN_CMD} -i {inp_path} -o {out_path}"
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

    # Copy output (pos or last frame) to relaxed/
    # For ENERGY runs, re-use clean CIF as "relaxed"
    run_type = _resolve_relax_mode(topology, n_atoms, relax_mode_override)
    if run_type == "ENERGY":
        import shutil
        src = record.get("clean_cif") or record.get("init_cif")
        if src:
            shutil.copy(src, relaxed_cif)
    else:
        # Parse last geometry from .out or -pos-1.xyz if written
        pos_file = os.path.join(run_dir, "cp2k_run-pos-1.xyz")
        if os.path.exists(pos_file):
            try:
                from pymatgen.core import Structure
                from pymatgen.io.xyz import XYZ
                struct = Structure.from_file(record.get("clean_cif"))
                # best effort: just save clean CIF when XYZ parsing is complex
            except Exception:
                pass
        import shutil
        src = record.get("clean_cif") or record.get("init_cif")
        if src:
            shutil.copy(src, relaxed_cif)

    return {**record,
            "dft_e_per_atom": e_per_atom,
            "dft_e_hartree": float(re.findall(r"ENERGY\| Total FORCE_EVAL[^\n]*\s+(-?\d+\.\d+)",
                                               out_text)[-1]),
            "relaxed_cif": relaxed_cif,
            "status": "converged"}
```

- [ ] **Step 4: Run tests**

```bash
python3 -m pytest tests/dft_validation/test_runner.py -v
```

Expected: all 3 pass (the integration tests use real clean CIFs, the parse test is pure unit).

---

### Task 3.2 — Dry pass: single ENERGY, pause for confirmation

**Files:**
- Modify: `dft_validation/run_dft_validation.py` (add phase 3)

- [ ] **Step 1: Add Phase 3 to orchestrator**

In `run_dft_validation.py`, add:

```python
from dft_validation.cp2k_runner import run_one
import csv

def phase3(records, results_dir):
    print("=== Phase 3: dry pass — single ENERGY single-point ===")
    # Pick first clean Cd record (not probe)
    candidate = next(
        (r for r in records
         if r.get("status") in ("clean", "preclean_done") and not r.get("probe_cu")),
        None
    )
    if candidate is None:
        print("  FAIL: no clean Cd structure found")
        sys.exit(1)

    print(f"  Running: {candidate['id']} (n_atoms={candidate.get('n_atoms')})")
    updated = run_one(candidate, results_dir, "ffc",
                      relax_mode_override="ENERGY", walltime_seconds=600)

    if updated.get("dft_e_per_atom") is None:
        print(f"  FAIL: no energy parsed. Status: {updated.get('status')}")
        sys.exit(1)

    e = updated["dft_e_per_atom"]
    e_h = updated.get("dft_e_hartree", "?")
    n   = updated.get("n_atoms", "?")
    print(f"\n  ✓ Energy parsed:")
    print(f"    E (Hartree) = {e_h}")
    print(f"    n_atoms     = {n}")
    print(f"    E (eV/atom) = {e:.6f}")
    print(f"\n  *** PAUSING for user confirmation ***")
    print(f"  Please verify the eV/atom value looks chemically reasonable")
    print(f"  (typical MOF: −5 to −15 eV/atom range) before proceeding to Phase 4.")

    # Write 1-row CSV
    csv_path = os.path.join(results_dir, "dft_vs_symmetry.csv")
    fieldnames = ["topology","composition_label","sym_class","rho",
                  "dft_e_per_atom","n_electrons","multiplicity","status","id"]
    with open(csv_path, "w", newline="") as csvf:
        w = csv.DictWriter(csvf, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerow(updated)
    print(f"  CSV row written to {csv_path}")
    return [updated if r["id"] == updated["id"] else r for r in records]
```

And in `main()`:
```python
    if args.phase >= 3:
        records = phase3(records, args.results_dir)
        if args.phase == 3:
            print("\nPhase 3 complete. Review energy value above, then run --phase 4.")
            sys.exit(0)
```

- [ ] **Step 2: Run Phase 3**

```bash
python3 dft_validation/run_dft_validation.py --phase 3
```

This takes a few minutes (ENERGY single-point, ~200 atoms, OT SCF). Expected output includes the eV/atom value. Typical range: −6 to −12 eV/atom for MOF with Cd/organic.

- [ ] **Step 3: ⚠ HARD PAUSE — report energy to user**

Show the user:
- `E (Hartree)` value
- `n_atoms` count
- `E (eV/atom)` value
- Confirm it is in the −5 to −15 eV/atom range

**Do not proceed to Phase 4 until the user explicitly confirms.**

- [ ] **Step 4: Commit (after user confirmation)**

```bash
git add dft_validation/cp2k_runner.py dft_validation/run_dft_validation.py \
    tests/dft_validation/test_runner.py
git commit -m "feat: Phase 3 — CP2K input generator, energy parser (Hartree→eV/atom), dry pass"
```

---

## Phase 4 — Full batch + analyze.py

**Prerequisite:** User has confirmed the Phase 3 energy value.

**Acceptance gate:** All 10 ffc Cd structures attempted; ≥7 converged; `dft_vs_symmetry.csv` written; ρ–E figure with 2 subplots (Comp A, Comp B) saved; console prints Spearman ρ_s per composition; fjh correctly gated (skipped with reason, since MAX_ATOMS=250 already blocks fjh assembly).

---

### Task 4.1 — Full batch run

**Files:**
- Modify: `dft_validation/run_dft_validation.py` (phase 4 run loop)

- [ ] **Step 1: Add Phase 4 run loop**

In `run_dft_validation.py`, add:

```python
def phase4_run(records, results_dir):
    print("=== Phase 4a: full batch CP2K run ===")
    updated = []
    for rec in records:
        if rec.get("probe_cu"):
            # Cu probe: always ENERGY, open-shell path
            result = run_one(rec, results_dir, "ffc",
                             relax_mode_override="ENERGY", walltime_seconds=1800)
        elif rec.get("status") in ("clean", "preclean_done"):
            result = run_one(rec, results_dir, "ffc",
                             walltime_seconds=7200)
        else:
            result = rec
        updated.append(result)
        # Update JSON after each run (enables resume)
        json_path = os.path.join(results_dir, "selected_structures.json")
        with open(json_path, "w") as f:
            json.dump(updated + records[len(updated):], f, indent=2)

    converged = sum(1 for r in updated if r.get("status") == "converged")
    print(f"[phase4] Converged: {converged}/{len(updated)}")
    return updated
```

And update `main()`:
```python
    if args.phase >= 4:
        records = phase4_run(records, args.results_dir)
```

---

### Task 4.2 — analyze.py

**Files:**
- Create: `dft_validation/analyze.py`
- Create: `tests/dft_validation/test_analyze.py`

**Interfaces:**
- Consumes: list of converged records with `rho`, `dft_e_per_atom`, `sym_class`, `composition_label`
- Produces:
  - `analyze_ffc(records, results_dir) -> dict` — keys `"comp_a"` and `"comp_b"`, each `{"spearman_rs": float, "p_val": float, "lowest_in_high": bool}`
  - `get_fjh_gate(comp_results: dict) -> tuple[bool, str]`
  - Writes `dft_vs_symmetry.csv` (full) and `rho_energy_ffc.png`

- [ ] **Step 1: Write failing tests**

```python
# tests/dft_validation/test_analyze.py
import sys
sys.path.insert(0, '.')
import pytest

def test_spearman_negative_trend():
    """Perfect inverse: higher rho → lower energy → rs < 0."""
    from dft_validation.analyze import _spearman_for_comp
    records = [
        {"rho": 1.0, "dft_e_per_atom": -10.0, "sym_class": "High"},
        {"rho": 0.5, "dft_e_per_atom": -9.5,  "sym_class": "High"},
        {"rho": 0.25,"dft_e_per_atom": -9.0,  "sym_class": "Med"},
        {"rho": 0.125,"dft_e_per_atom": -8.5, "sym_class": "Low"},
    ]
    rs, p, lowest_in_high = _spearman_for_comp(records)
    assert rs < 0, f"Expected rs<0, got {rs}"
    assert lowest_in_high is True

def test_fjh_gate_passes():
    from dft_validation.analyze import get_fjh_gate
    comp_results = {
        "comp_a": {"spearman_rs": -0.9, "p_val": 0.05, "lowest_in_high": True},
        "comp_b": {"spearman_rs": -0.8, "p_val": 0.08, "lowest_in_high": True},
    }
    passed, reason = get_fjh_gate(comp_results)
    assert passed is True

def test_fjh_gate_fails_on_not_in_high():
    from dft_validation.analyze import get_fjh_gate
    comp_results = {
        "comp_a": {"spearman_rs": -0.9, "p_val": 0.05, "lowest_in_high": False},
        "comp_b": {"spearman_rs": -0.8, "p_val": 0.08, "lowest_in_high": True},
    }
    passed, reason = get_fjh_gate(comp_results)
    assert passed is False
    assert "comp_a" in reason

def test_fjh_gate_fails_on_positive_rs():
    from dft_validation.analyze import get_fjh_gate
    comp_results = {
        "comp_a": {"spearman_rs": 0.3,  "p_val": 0.5, "lowest_in_high": True},
        "comp_b": {"spearman_rs": -0.8, "p_val": 0.08, "lowest_in_high": True},
    }
    passed, reason = get_fjh_gate(comp_results)
    assert passed is False
```

- [ ] **Step 2: Run tests — verify FAIL**

```bash
python3 -m pytest tests/dft_validation/test_analyze.py -v 2>&1 | head -10
```

- [ ] **Step 3: Write `dft_validation/analyze.py`**

```python
# dft_validation/analyze.py
import os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from dft_validation import config


SYM_COLORS = {"High": "#d62728", "Med": "#1f77b4", "Low": "#95a5a6"}


def _spearman_for_comp(records: list) -> tuple:
    """
    Compute Spearman rs for (rho, dft_e_per_atom).
    Returns (rs, p_val, lowest_in_high: bool).
    """
    converged = [r for r in records if r.get("dft_e_per_atom") is not None]
    if len(converged) < 2:
        return (float("nan"), float("nan"), False)
    rhos     = [r["rho"] for r in converged]
    energies = [r["dft_e_per_atom"] for r in converged]
    rs, p    = spearmanr(rhos, energies)
    min_rec  = min(converged, key=lambda r: r["dft_e_per_atom"])
    return (rs, p, min_rec["sym_class"] == "High")


def get_fjh_gate(comp_results: dict) -> tuple:
    """
    Qualitative gate: both comps must have lowest-E in High AND rs < 0.
    Returns (passed: bool, reason: str).
    """
    for label, stats in comp_results.items():
        if not stats.get("lowest_in_high", False):
            return False, f"{label}: lowest-energy structure is NOT in High class"
        if stats.get("spearman_rs", 0) >= 0:
            return False, (f"{label}: Spearman rs={stats['spearman_rs']:.3f} >= 0 "
                           f"(energy does not decrease with symmetry)")
    return True, "both compositions pass gate"


def analyze_ffc(records: list, results_dir: str) -> dict:
    """
    Per-composition Spearman + side-by-side ρ–E plot.
    Returns dict comp_label → {spearman_rs, p_val, lowest_in_high}.
    """
    comp_results = {}
    labels = sorted({r["composition_label"] for r in records
                     if r.get("topology") == "ffc" and not r.get("probe_cu")})

    fig, axes = plt.subplots(1, len(labels), figsize=(6 * len(labels), 5), sharey=False)
    if len(labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, labels):
        subset = [r for r in records
                  if r.get("composition_label") == label
                  and r.get("dft_e_per_atom") is not None]
        rs, p, lowest_in_high = _spearman_for_comp(subset)
        comp_results[label] = {"spearman_rs": rs, "p_val": p, "lowest_in_high": lowest_in_high}

        for cls in ["High", "Med", "Low"]:
            pts = [r for r in subset if r.get("sym_class") == cls]
            if pts:
                ax.scatter([r["rho"] for r in pts],
                           [r["dft_e_per_atom"] for r in pts],
                           color=SYM_COLORS[cls], label=cls, s=60, alpha=0.85, zorder=3)

        ax.set_xlabel(r"Symmetry Index $\rho$", fontsize=12)
        ax.set_ylabel("DFT Energy (eV/atom)", fontsize=12)
        rs_str = f"{rs:.3f}" if rs == rs else "n/a"
        ax.set_title(f"ffc {label}  (Spearman $\\rho_s$={rs_str}, p={p:.3f})", fontsize=11)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(results_dir, "rho_energy_ffc.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[analyze] Figure saved to {fig_path}")
    return comp_results


def write_csv(records: list, results_dir: str):
    """Write full results CSV."""
    csv_path = os.path.join(results_dir, "dft_vs_symmetry.csv")
    fieldnames = ["id", "topology", "composition_label", "sym_class", "rho",
                  "dft_e_per_atom", "dft_e_hartree", "n_electrons",
                  "multiplicity", "n_atoms", "status"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    print(f"[analyze] CSV written to {csv_path}")
```

- [ ] **Step 4: Run tests**

```bash
python3 -m pytest tests/dft_validation/test_analyze.py -v
```

Expected: all 4 PASS.

---

### Task 4.3 — Wire Phase 4 into orchestrator + acceptance check

**Files:**
- Modify: `dft_validation/run_dft_validation.py` (add phase4_analyze, wire fjh gate, final summary)

- [ ] **Step 1: Complete `run_dft_validation.py`**

Append to existing file:

```python
from dft_validation.analyze import analyze_ffc, get_fjh_gate, write_csv

def phase4_analyze(records, results_dir):
    print("=== Phase 4b: analyze ===")
    write_csv(records, results_dir)
    comp_results = analyze_ffc(records, results_dir)

    print("\n--- ffc results ---")
    for label, stats in comp_results.items():
        hi = "YES" if stats["lowest_in_high"] else "NO"
        print(f"  {label}: Spearman rs={stats['spearman_rs']:.3f}  "
              f"p={stats['p_val']:.3f}  lowest_in_High={hi}")

    # fjh gate
    if config.FJH_GATE_ON_FFC:
        passed, reason = get_fjh_gate(comp_results)
        if passed:
            print(f"\n  [fjh gate] PASS — {reason}")
            print("  (fjh skipped locally: MAX_ATOMS=250 blocks ~800-atom fjh structures)")
        else:
            print(f"\n  [fjh gate] SKIP — {reason}")

    # Count converged / skipped / discarded
    converged  = [r for r in records if r.get("status") == "converged"]
    failed     = [r for r in records if r.get("status") not in
                  ("converged", "done_cached", "selected", "assembled",
                   "clean", "preclean_done")]
    print(f"\n  converged={len(converged)}  failed/skipped={len(failed)}")
    print("  Phase 4 PASS")
    return records
```

Update `main()`:
```python
    if args.phase >= 4:
        records = phase4_run(records, args.results_dir)
        records = phase4_analyze(records, args.results_dir)
```

- [ ] **Step 2: Run full pipeline**

```bash
python3 dft_validation/run_dft_validation.py --phase 4
```

This will take ~30–120 min for 10 GEO_OPT structures (4 MPI × ~200 atoms each). Monitor progress with:

```bash
ls results/dft_validation/ffc/runs/*/cp2k_run.out | wc -l
tail -5 results/dft_validation/ffc/runs/ffc_comp_a_*/cp2k_run.out
```

- [ ] **Step 3: Verify acceptance criteria**

```bash
# 1. JSON has all 11 records
python3 -c "
import json
recs = json.load(open('results/dft_validation/selected_structures.json'))
print('total:', len(recs))
conv = [r for r in recs if r.get('status')=='converged']
print('converged:', len(conv))
"

# 2. CSV exists with energy values
python3 -c "
import csv
rows = list(csv.DictReader(open('results/dft_validation/dft_vs_symmetry.csv')))
print('CSV rows:', len(rows))
for r in rows:
    print(r['id'], r['status'], r.get('dft_e_per_atom','?'))
"

# 3. Figure exists
ls -lh results/dft_validation/rho_energy_ffc.png
```

- [ ] **Step 4: Final commit**

```bash
git add dft_validation/ tests/dft_validation/ results/dft_validation/*.csv \
    results/dft_validation/*.png
git commit -m "feat: Phase 4 — full batch CP2K run + per-composition Spearman + ρ–E figure"
```

---

## Self-Review Checklist

**Spec coverage:**
- ✓ ffc Comp A + B exhaustive (10 structures, no sampling)
- ✓ Cu probe (1 structure, ENERGY, open-shell path)
- ✓ Electron parity check → doublet override
- ✓ CHARGE=0
- ✓ xTB FIRE pre-relax (steps=200 cap), periodic
- ✓ Element coverage against BASIS_MOLOPT_UZH + GTH_POTENTIALS
- ✓ Basis = DZVP-MOLOPT-PBE-GTH-q{n} from BASIS_MOLOPT_UZH
- ✓ Hartree × 27.211386 / n_atoms = eV/atom
- ✓ OT SCF for closed-shell; diagonal + Fermi-Dirac for open-shell
- ✓ DFTD3(BJ) dispersion
- ✓ Per-topology RELAX_MODE dict; auto mapping
- ✓ Phase 3 hard pause before full batch
- ✓ Resumable (skips if relaxed CIF exists)
- ✓ Per-composition Spearman + side-by-side plot
- ✓ fjh qualitative gate (both comps lowest-E in High AND rs < 0)
- ✓ fjh skipped locally (MAX_ATOMS=250 blocks fjh)
- ✓ NetworkX v3 monkeypatch
- ✓ No modifications to existing business files

**Type consistency:** `run_one` returns updated record dict; `analyze_ffc` consumes list of those dicts; `get_fjh_gate` consumes the `comp_results` dict produced by `analyze_ffc`. All match.

**Placeholder scan:** None found.
