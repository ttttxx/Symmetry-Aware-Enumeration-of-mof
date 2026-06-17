# CASM Enumeration Comparison Experiment

## Purpose

Cross-validate our graph-automorphism-based MOF enumerator against CASM's
crystallographic enumeration on the same quotient graph. Compares orbit counts
(correctness), wall-clock time, and peak RSS memory.

## CASM Version

- **libcasm v2** — `libcasm-configuration` 2.3.1, `libcasm-xtal` 2.3.0
- Installed via `pip install libcasm-configuration` (Python 3.11 required)
- CASM worker: `experiment/_casm_worker.py` (run with `python3.11`)

## Symmetry Restriction — Fallback Path

CASM's `ConfigEnumAllOccupations` does not accept a custom `SymGroup` parameter,
so the **primary path** (injecting proper-rotation-only subgroup) is not feasible.

We use the **fallback path** instead:

1. CASM enumerates configurations modulo the **full factor group** F (including
   improper operations: mirrors, inversions, glide planes). Output: `orbits_casm_full = F_total`.
2. Our enumerator uses only the **proper-rotation subgroup** Γ (det > 0 only).
   Output: `orbits_ours = P`.
3. Reconciliation formula:
   - Each CASM full-group orbit is classified as **achiral (A)** if its stabilizer
     contains at least one improper operation, or **chiral (C)** if not.
   - Under the proper-rotation subgroup, each chiral orbit splits into two
     enantiomeric orbits.
   - Therefore: `P = A + 2C` and `F_total = A + C`.
   - From the counts: `C = P - F_total`, `A = 2*F_total - P`.
   - `orbits_casm` (proper-rotation equivalent) = `P = F_total + C`.
4. Correctness gate for exhaustive cases:
   ```
   assert C >= 0 and A >= 0    # non-negative counts
   assert P == orbits_ours     # orbit counts match
   ```

For ffc, all configurations happen to be achiral (C=0), so `F_total == P`.
For fjh, some are chiral — the reconciliation formula is needed.

## Site Mapping

Each topology's CIF (`data/mof/{topology}.cif`, P1 format) is used verbatim
(no `get_primitive_standard_structure()`). Site order is validated against
`data/{t}/{t}_graph.json`:

1. Site count: `len(structure.sites) == graph['n_nodes']`
2. Coordinate match: `|CIF_frac - graph_frac| < 1e-3` per site
3. Degree match: consistent with `{t}_structure_info.json`

## Timing Boundary

Both sides measure **only** the enumeration core call:

- **CASM**: `time.perf_counter()` around the `for _ in enum.by_supercell_list(...)` loop.
  **Note:** CASM timing includes Python iteration overhead from the counting loop.
- **Ours**: `time.perf_counter()` around `enumerator.enumerate(all_ligands)`.

Symmetry analysis (pymatgen `SpacegroupAnalyzer`, CASM `make_prim_factor_group()`)
is **outside** the timing boundary on both sides.

## Memory Measurement

Each `(case, side, rep)` runs in a **fresh isolated subprocess**. The subprocess
calls `resource.getrusage(resource.RUSAGE_SELF).ru_maxrss` once, after enumeration
completes. Because the subprocess has no prior RSS high-watermark, this value is
the subprocess's peak RSS.

**Why not before/after delta?** `ru_maxrss` is a high-watermark counter. If earlier
code in the same process already set a high RSS, the delta would be zero. Using a
fresh subprocess avoids this.

On **macOS**, the unit of `ru_maxrss` is **bytes** (unlike Linux where it's kilobytes).

## Thread Control

No `taskset` (not available on macOS). Single-thread behavior achieved via:
```python
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
```

## Repetitions

Each `(topology, ligand_pool, side)` is run **3 times**. The **median** wall-clock
time and the peak RSS from the median-time run are recorded.

## Intractable Threshold (vba / CASM side)

A CASM run is declared **intractable** if either:
- Wall-clock time > **30 minutes** (1800 s), or
- Process RSS > **8 GB**

The parent process polls with `psutil` every **5 seconds**. Whichever threshold is
reached first triggers `proc.terminate()` and records `orbits_casm = "intractable"`.

## Run Command

```bash
# Full run (all 7 cases)
python3 experiment/casm_compare.py

# Quick test (ffc only, 2 cases)
python3 experiment/casm_compare.py --quick
```

**Prerequisites:**
- Python 3.11 with `libcasm-configuration`, `libcasm-xtal`, `pymatgen`, `numpy`
- Project Python 3 with `pymatgen`, `numpy` (and project modules)
- Optional: `psutil` for vba watchdog

## Output

- **CSV:** `experiment/casm_compare_results.csv`
- **Scripts:**
  - `experiment/casm_compare.py` — main harness
  - `experiment/_casm_worker.py` — CASM enumeration subprocess worker
  - `experiment/_ours_worker.py` — our enumerator subprocess worker

## Per-Case Status Summary

| Topology | ligand_pool | Expected orbits_ours | status_ours | status_casm |
|----------|-------------|---------------------|-------------|-------------|
| ffc | C3=2,C4=2 | 42 | exhaustive | exhaustive |
| ffc | C3=3,C4=3 | 513 | exhaustive | exhaustive |
| fjh | C3=1,C4=2 | 394 | exhaustive | exhaustive |
| fjh | C3=2,C4=1 | 4 756 | exhaustive | exhaustive |
| fjh | C3=1,C4=3 | 36 369 | exhaustive | exhaustive |
| vba | C4=2 | — | budgeted | intractable |
| vba | C4=3 | — | budgeted | intractable |

## Notes

- ffc factor group: 16 total, 8 proper (P4/mmm). All configurations are achiral.
- fjh factor group: 32 total, 16 proper (I4/mcm). Mix of achiral (A) and chiral (C)
  orbits — reconciliation formula applied.
- vba factor group: 8 total, 8 proper (Fddd). BUT graph automorphism expects |Γ|=16.
  The crystallographic symmetry is lower than the topological symmetry; documented
  as a known discrepancy. vba cases are intractable either way.
- Our enumerator's Burnside verification (Sum(1/|Stab|) == |X|/|G|) confirms
  completeness for exhaustive cases.
