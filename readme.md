# Symmetry-Aware Enumeration and Reconstruction for MOF Candidate Generation

This repository contains the research code accompanying the manuscript:

**From Topological Design to Physically Instantiable MOF Candidates: Symmetry-Aware Enumeration and Reconstruction**

The codebase implements a deterministic workflow for:

1. loading a periodic MOF topology as a quotient-graph-style graph representation,
2. extracting a proper-rotation symmetry group from that topology,
3. enumerating symmetry-distinct building-block assignments,
4. verifying enumeration completeness with an orbit-stabilizer / Burnside identity when the run is exhaustive,
5. reconstructing accepted assignments into periodic crystal models, and
6. optionally carrying representative candidates into a screening-stage DFT pipeline.

This repository is a research-code release aligned with the manuscript. The main purpose of this README is to explain how to set up the environment and run the core workflow.

## Quick Start

If you only want the shortest path to a working run, use the following steps.

### 1. Create and activate a virtual environment

Using a virtual environment is strongly recommended. The code depends on scientific Python packages such as `pymatgen`, `ase`, and `networkx`, so installing them into an isolated environment is the safest option.

macOS / Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Windows:

```bash
python -m venv .venv
.venv\Scripts\activate
```

If the environment is active, your shell prompt usually shows something like `(.venv)`.

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

Optional import check:

```bash
python test_import.py
```

### 3. One-command run

The main entry point of this submission version is:

```bash
python main.py
```

This command runs the built-in example configured at the bottom of `main.py` and performs:

1. topology loading and enrichment,
2. symmetry extraction,
3. symmetry-distinct enumeration,
4. constructive verification,
5. reconstruction of representative structures into CIF files.

Generated files are written under `results/`.

### 4. Exhaustive run where feasible

If you want to disable the default enumeration budgets:

```bash
python main.py --disable-budget
```

### 5. Optional benchmark JSON

```bash
python main.py --benchmark-json results/performance/demo.json
```

## What This Repository Covers

### 1. Topology to finite design space

Prepared topology graphs are stored under `data/<topology>/<topology>_graph.json`. The main pipeline enriches these graphs with periodic-topology metadata needed for reconstruction and cycle-based bookkeeping.

### 2. Symmetry extraction

`symmetry.py` uses `pymatgen` symmetry operations and retains only proper rotations, because the manuscript treats symmetry reduction as reconstruction-compatible rather than as a full crystallographic factor-group quotient.

### 3. Symmetry-distinct enumeration

`enumerator.py` organizes the search by subgroup-induced orbit partitions, applies degree-aware and count-aware pruning, and stores canonical representatives in a trie. In exhaustive mode, the output can be checked exactly against the counting identity.

### 4. Constructive self-verification

`verification.py` computes the constructive sum

```text
sum_s 1 / |Stab(s)|
```

and compares it with the theoretical counting value

```text
|Fix(e)| / |G|
```

to determine whether the enumerated representatives close the full orbit space.

### 5. Reconstruction into crystal models

`assembler.py` maps accepted assignments back to periodic framework models by:

- adaptive lattice scaling,
- node placement,
- Kabsch-based rigid-body alignment,
- linker placement,
- periodicity-aware edge adjustment, and
- CIF export plus structural validation.

### 6. Reconstruction and DFT follow-up

The repository also includes:

- `experiment/geometric_reconstruction_summary.py` for representative reconstruction audits across multiple topologies,
- `dft_validation/` for selection, sanitization, CP2K execution, and analysis of screening-stage DFT examples.

## Repository Layout

```text
mof/
├── main.py                              # End-to-end enumeration + verification + CIF generation
├── symmetry.py                          # Proper-rotation symmetry extraction
├── enumerator.py                        # Symmetry-distinct enumeration with pruning and budgets
├── verification.py                      # Burnside / orbit-stabilizer self-verification
├── assembler.py                         # Rigid-body reconstruction into periodic CIFs
├── data/                                # Prepared topology graph inputs and source CIFs
├── nodes/                               # Node building-block CIF library
├── edges/                               # Edge / linker CIF library
├── experiment/                          # Benchmarking and manuscript-facing audit scripts
├── dft_validation/                      # Screening-stage DFT workflow
├── tests/                               # Tests for selected components
├── external_lib/                        # Reused periodic-geometry utilities
├── results/                             # Generated outputs (enumeration, CIFs, DFT, audits)
└── paper/                               # Manuscript sources and figures
```

## Environment Setup

The code was developed in Python 3.8+.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Optional import smoke test:

```bash
python test_import.py
```

Core dependencies for the manuscript-aligned CLI workflow include `numpy`, `pymatgen`, `networkx`, `ase`, and `matplotlib`. `chgnet` is only needed for optional parts of the validation workflow.

For this manuscript-aligned submission version, the core workflow does not depend on `app.py`. The command-line pipeline centered on `main.py` is the primary supported entry point.

## Reproducible Workflows

### 1. End-to-end enumeration and reconstruction

`main.py` performs the following stages:

1. load a prepared topology graph,
2. enrich periodic-topology metadata,
3. compute the automorphism group,
4. enumerate symmetry-distinct assignments,
5. verify completeness when applicable,
6. reconstruct representative assignments into CIFs,
7. write JSON summaries to `results/`.

### 2. Reconstruction audit used for manuscript-facing summaries

```bash
python experiment/geometric_reconstruction_summary.py
```

Useful variants:

```bash
python experiment/geometric_reconstruction_summary.py --topology ffc --topology fjh
python experiment/geometric_reconstruction_summary.py --per-class 3
```

This workflow samples representative structures by symmetry class and reports assembly success, sanitization status, topology recovery, and geometric closure metrics.

### 3. Screening-stage DFT workflow

The DFT pipeline lives in `dft_validation/` and is driven by `dft_validation/run_dft_validation.py`.

Typical commands are:

```bash
python -u dft_validation/run_dft_validation.py --phase 1
python -u dft_validation/run_dft_validation.py --phase 2
python -u dft_validation/run_dft_validation.py --phase 3
python -u dft_validation/run_dft_validation.py --phase 4 --summary
python -u dft_validation/run_dft_validation.py --resume --harvest-only --summary
python -u dft_validation/run_dft_validation.py --resume --analyze-only --summary
```

Please review `dft_validation/config.py` before running this pipeline. It contains topology targets, CP2K settings, memory limits, and execution-mode switches.

## Main Outputs

Depending on the workflow, the repository produces:

- enumeration JSON files in `results/`,
- reconstructed CIF directories in `results/*_cifs/`,
- reconstruction audit CSV/JSON summaries in `results/reconstruction_summary*` and related folders,
- DFT selections, parsed energies, and plots in `results/dft_validation/`.

## Important Interpretation Boundaries

This repository should be read together with the manuscript's stated limits.

- The reconstruction stage is a rigid-body workflow. It is appropriate for rigid or semi-rigid building units, but it is not a conformational search engine for highly flexible linkers.
- The DFT workflow is a screening-stage example based on fixed-composition comparisons and single-point calculations on as-assembled structures. It is not a final thermodynamic ranking pipeline across chemistries.
- Larger search spaces may be explored with explicit runtime or result budgets. In those cases, the output should be interpreted as partial coverage unless the verification identity closes exactly.
- The five-topology reconstruction gallery and follow-up audits are representative validation steps, not a universal success-rate claim for all MOF topologies.

## Paper-to-Code Map

For readers matching the manuscript to the implementation:

- `Candidate Generation and Reconstruction` -> `main.py`, `symmetry.py`, `enumerator.py`, `assembler.py`
- `Algebraic Correctness and External Cross-Check` -> `verification.py`, `experiment/casm_compare.py`
- `Geometric Instantiation Across Representative Topologies` -> `experiment/geometric_reconstruction_summary.py`
- `Illustrative Screening-Stage DFT Example` -> `dft_validation/`

## Notes on Repository State

This repository currently contains research artifacts in addition to source code, including intermediate results, manuscript files, and validation outputs. For manuscript review, the most relevant source files are the Python modules listed above rather than every file under `results/` or `paper/`.
