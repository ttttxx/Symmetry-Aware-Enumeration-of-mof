# CASM Enumeration Comparison Experiment — Design Spec

**Date:** 2026-06-10
**Status:** Approved (with corrections applied)

---

## Goal

Run CASM's algebraic site-occupation enumeration on the same quotient graph as our enumerator,
producing cross-validated orbit counts plus wall-clock time and peak RSS for a direct performance
comparison. The comparison covers five exhaustive benchmark cases (ffc × 2, fjh × 3) and two
budget-truncated non-exhaustive cases (vba × 2).

---

## Scientific Constraints (must not be violated)

1. **Proper-rotation-only symmetry group.** Our Γ contains only proper rotations
   (Cartesian rotation matrix det = +1). Mirrors, inversions, and glide planes are excluded;
   enantiomeric configurations are counted as distinct. CASM's factor group must be filtered to
   the proper-rotation subgroup before any enumeration begins.

2. **Degree-constrained occupation space.** Each site's allowed species = the ligand pool for
   that site's degree. The enumeration space is `∏_d |C_d|^{n_d}` (product over degree classes).
   Supercell volume = 1 throughout; no supercell expansion.

3. **Correctness gate.** Exhaustive cases must match reference orbit counts exactly:

   | Topology | ligand_pool | Expected orbits |
   |----------|-------------|-----------------|
   | ffc | C3=2,C4=2 | 42 |
   | ffc | C3=3,C4=3 | 513 |
   | fjh | C3=1,C4=2 | 394 |
   | fjh | C3=2,C4=1 | 4 756 |
   | fjh | C3=1,C4=3 | 36 369 |

   Any mismatch → stop, report, do not proceed.

---

## Allowed Code Changes

| File | Change |
|------|--------|
| `enumerator.py` | Add `--benchmark-json` output port only (non-invasive) |
| `main.py` | Add `argparse` flag + timing/memory wrapper calling the same `enumerator.enumerate()` |
| `experiment/casm_compare.py` | New script (CASM prim build, symmetry filter, enumeration, harness) |
| `experiment/casm_compare_README.md` | New documentation |

No changes to enumeration logic, pruning, deduplication, or any code that affects orbit counts.

---

## Step 1 — Non-invasive Benchmark Output Port

### CLI flag

```
python main.py --benchmark-json results/bench_ffc_c3x2c4x2.json
```

Default off; absent flag → behaviour identical to current.

### Timing boundary

`time.perf_counter()` wraps **only** `enumerator.enumerate(all_ligands)` — excludes
pymatgen symmetry extraction and Kabsch reconstruction. This matches CASM's enumeration
boundary defined in Step 4.

### Memory boundary

`resource.getrusage(resource.RUSAGE_SELF).ru_maxrss` is read immediately before and after
`enumerate()`. On **macOS** the unit is **bytes** (not kilobytes). The delta from the
pre-measurement baseline is stored; if the OS did not update RSS during the call the
post-call value is used directly. `tracemalloc` is **not** used.

### JSON record schema

```json
{
  "topology": "ffc",
  "ligand_pool": "C3=2,C4=2",
  "orbits": 42,
  "time_s": 0.031,
  "peak_mem_mb": 12.4,
  "exhaustive": true,
  "budget_triggered": false
}
```

`ligand_pool` encoding: `"C{d}={k}"` pairs sorted by degree, comma-separated
(e.g. `"C3=1,C4=3"`). This ensures alignment with CASM output rows.

`budget_triggered` is `True` when `task_time_limit` cut off any subgroup task **or**
`IDENTITY_NEW_LIMIT` / `IDENTITY_ORBIT_LIMIT` was hit.

### Regression gate

After adding the port, run ffc and fjh exhaustive cases with the flag **off** and confirm
orbit counts remain 42 / 513 / 394 / 4 756 / 36 369 before trusting any timing numbers.

---

## Step 2 — CASM Prim Construction

### Source CIF

Use `data/mof/{topology}.cif` verbatim — **no `get_primitive_standard_structure()`**, no
re-standardisation, no transformation of site order. The cell is stored as P1 (symmetry
stripped); space-group symmetry is re-derived in Step 3.

### Site mapping validation (required, explicit)

Before building the CASM prim, run three assertions against `data/{t}/{t}_graph.json`:

1. **Site count**: `len(structure.sites) == graph_data['n_nodes']`
2. **Coordinate match**: fractional coordinates of each CIF site agree with the node's
   `frac` field in the graph JSON (within tolerance 1e-3).
3. **Degree match**: the site's element symbol maps to the expected degree class.
   Mapping is read from `data/{t}/{t}_structure_info.json` which stores `species` and the
   graph's `node_degrees` list.

Failure of any assertion → raise `ValueError` with a descriptive message; do not proceed.

### Species assignment

For a given benchmark case `(topology, |C3|, |C4|)`:

- degree-3 sites → allowed species `["A", "B", ...]` (first `|C3|` uppercase letters)
- degree-4 sites → allowed species `["X", "Y", ...]` (last `|C4|` uppercase letters, starting from X)

The placeholder symbols carry no physical meaning; they exist only for combinatorial counting.

### libcasm prim construction

```python
from libcasm.crystallography import BasicStructure, Site, Coordinate

prim = BasicStructure(
    lattice=Lattice(structure.lattice.matrix),
    title="ffc_compare"
)
for i, site in enumerate(structure.sites):
    degree = node_degrees[i]
    species = species_for_degree(degree, c3_count, c4_count)
    prim.push_back(Site(Coordinate(site.frac_coords, prim.lattice), species))
```

---

## Step 3 — Proper-Rotation Subgroup

### Primary path (timebox 30 min)

1. Load the prim's factor group:
   ```python
   fg = prim.factor_group()  # list of SymOp
   ```
2. For each `op` extract the Cartesian rotation matrix `R`:
   ```python
   R = op.matrix  # 3×3 float64
   ```
3. Filter:
   ```python
   proper_ops = [op for op in fg if np.linalg.det(op.matrix) > 0]
   ```
   Screw axes (proper rotation + fractional translation) are **kept** because
   `det(R) > 0` regardless of translation. Mirrors, glide planes, and inversion are dropped.

4. **Mandatory assertion** (correctness gate):
   ```python
   expected = {"ffc": 8, "fjh": 16, "vba": 16}
   assert len(proper_ops) == expected[topology], (
       f"Proper-rotation subgroup size {len(proper_ops)} ≠ expected {expected[topology]} "
       f"for {topology}. Aborting."
   )
   ```
   If this fails → print a diagnostic table (all ops + their det values) and `sys.exit(1)`.

5. Construct CASM SymGroup from `proper_ops` and pass to `ConfigEnumAllOccupations`:
   ```python
   sym_group = SymGroup(proper_ops)
   enumerator = ConfigEnumAllOccupations(supercell, sym_group=sym_group)
   ```

### Fallback path (if libcasm API does not accept custom SymGroup)

If `ConfigEnumAllOccupations` does not accept `sym_group=`, activate fallback:

1. Run with the full factor group → `orbits_casm_full`
2. For each orbit representative `r`, compute its stabilizer under the full factor group.
3. Classify each orbit:
   - **A (achiral)**: stabilizer contains at least one improper op (det < 0). Stays as one
     orbit under Γ_proper.
   - **C (chiral)**: stabilizer contains only proper ops. Splits into two enantiomeric orbits
     under Γ_proper.
4. Notation: `F_total = A + C` (= `orbits_casm_full`)
5. Predicted proper-rotation count: `P = A + 2C`
6. Assertion:
   ```python
   assert P == F_total + C == orbits_ours
   ```
   If this fails → report A, C, P, orbits_ours and abort.
7. CSV gains extra column `orbits_casm_full`; README documents the reconciliation formula.

---

## Step 4 — Measurement Harness

### Thread control

Set before any numeric computation:
```python
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
```
`taskset` is **not** available on macOS; the above env-var approach achieves
single-BLAS-thread execution.

### Timing boundary

Both sides measure **only** the enumeration core call:

- **CASM**: `perf_counter()` around `list(ConfigEnumAllOccupations(...))` (forces full iteration)
- **Ours**: `perf_counter()` around `enumerator.enumerate(all_ligands)` (same as Step 1 port)

Symmetry analysis (pymatgen `SpacegroupAnalyzer`, CASM `prim.factor_group()`) is **outside**
the timing boundary on both sides.

### Memory boundary

Same function call boundaries as timing. On macOS:

```python
import resource
rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # bytes on macOS
[enumeration call]
rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss   # bytes on macOS
peak_mem_mb = (rss_after - rss_before) / 1_048_576
```

If `peak_mem_mb < 0` (RSS only increases), use `rss_after / 1_048_576` as absolute RSS.

### Repetitions

Each (topology, ligand_pool, side) combination is run **3 times**; the **median** wall-clock
time is recorded. Peak RSS is taken from the run with the median time.

### Exhaustive vs. budgeted configuration

- Exhaustive cases: `task_time_limit=None`, `IDENTITY_ORBIT_LIMIT=inf`, `IDENTITY_NEW_LIMIT=inf`
  on our side (disable all budget triggers). CASM has no budget concept.
- vba cases: our side runs with default budgets. CASM side runs until intractable threshold.

### Intractable threshold (vba / CASM side)

A CASM run is declared **intractable** if either:
- Wall-clock time > **30 minutes**, or
- Process RSS > **8 GB**

Whichever is reached first → kill the process, record `orbits_casm = "intractable"`,
`time_casm_s = "intractable"`, `mem_casm_mb = "intractable"`. The threshold values
are documented in `experiment/casm_compare_README.md`.

---

## Step 5 — Outputs

### Benchmark cases

| Row | topology | ligand_pool | status_ours | status_casm |
|-----|----------|-------------|-------------|-------------|
| 1 | ffc | C3=2,C4=2 | exhaustive | exhaustive |
| 2 | ffc | C3=3,C4=3 | exhaustive | exhaustive |
| 3 | fjh | C3=1,C4=2 | exhaustive | exhaustive |
| 4 | fjh | C3=2,C4=1 | exhaustive | exhaustive |
| 5 | fjh | C3=1,C4=3 | exhaustive | exhaustive |
| 6 | vba | C4=2 | budgeted | intractable |
| 7 | vba | C4=3 | budgeted | intractable |

### CSV schema (`experiment/casm_compare_results.csv`)

```
topology, ligand_pool,
orbits_ours, orbits_casm,
time_ours_s, time_casm_s,
mem_ours_mb, mem_casm_mb,
status_ours, status_casm
```

vba rows: `orbits_ours`, `time_ours_s`, `mem_ours_mb` contain the actual budgeted values
(not `"-"`). `orbits_casm`, `time_casm_s`, `mem_casm_mb` contain `"intractable"`.

### Script (`experiment/casm_compare.py`)

Single Python script. Sections in order:
1. Imports + env-var thread limits
2. `build_prim(topology, c3_count, c4_count)` — CIF load + site validation + prim construction
3. `get_proper_symgroup(prim, topology)` — factor group + det filter + assertion
4. `run_casm_case(topology, c3_count, c4_count, reps=3)` — measurement, returns dict
5. `run_ours_case(topology, ligands_dict, reps=3)` — calls `--benchmark-json` subprocess
6. `main()` — runs all 7 cases, writes CSV
7. `if __name__ == "__main__": main()`

### README (`experiment/casm_compare_README.md`)

Must document:
- libcasm version installed
- Symmetry restriction implementation (primary path or fallback, and which was used)
- Timing boundary definition (both sides)
- Memory boundary definition and macOS-specific unit note
- Thread-limit approach (env vars, no taskset)
- Intractable threshold values (30 min / 8 GB) and rationale
- Run command: `python experiment/casm_compare.py`
- Per-case status summary table

---

## Spec Self-Review

- **Placeholders**: All `「填:…」` items from task.txt resolved:
  - CIF paths: `data/mof/{ffc,fjh,vba}.cif` ✓
  - Degree→species mapping: degree-3 → A,B,…; degree-4 → X,Y,… ✓
  - Thread limiting: OMP/OPENBLAS/MKL env vars ✓
  - Repetition count N: 3 ✓
  - Intractable threshold: 30 min or 8 GB ✓
- **Internal consistency**: symmetry group sizes (ffc=8, fjh=16, vba=16) match
  `performance_experiment.py` line 118 comment.
- **Scope**: only `experiment/` additions + Step 1 port in `enumerator.py`/`main.py`.
- **Ambiguity**: fallback P=A+2C formula unambiguous; vba status_ours vs status_casm separated.
