# CASM / hybrid subgroup DFS scale comparison (2026-11-16)

Run from the project root:

```sh
.venv/bin/python tools/build_enumeration_core.py
.venv/bin/python -m pytest -q tests/test_casm_scale_benchmark.py tests/test_native_enumerator.py tests/test_enumerator_storage.py tests/test_counter_enumerator.py tests/test_subgroup_task_dedup.py tests/test_identity_budget.py
.venv/bin/python experiment/casm_scale_compare.py
.venv/bin/python validation/analyze_casm_scale_20261116.py
```

The parent uses the project environment with psutil. Both measurement workers
use **the same Python 3.11 interpreter**, which must already provide CASM,
NumPy, SymPy, pymatgen and project imports. Native ctypes does not depend on a
particular Python extension ABI. Override the worker interpreter with `--python`.

The 15 predefined cases are in `casm_scale_common.py`. Label pools are sorted
prefixes of existing `nodes/3c_*.cif` and `nodes/4c_*.cif` filenames. This is an
assignment-enumeration benchmark; it does not construct or validate new models
from all of those templates.

Four separately timed sides:

- `ours`: required C++ DFS, Python subgroup scheduling, full byte-key index,
  result dictionaries, stabilizers, orbit sizes, sorting and label decoding.
- `casm_records`: CASM full-candidate generator filtered by the exactly matched
  proper subgroup, with the same output fields, stabilizer calculation and sorting.
- `casm_count`: same proper subgroup but counting only, without retaining records.
- `casm_full`: the original CASM public API counting under its full factor group;
  its group differs, so it is a supplementary reference, not the matched speedup.

Each side runs three times in independent sequential processes. Graph/group
construction and dependency imports occur before timing; our subgroup-task
preparation remains inside timing. Both worker paths set OMP/OpenBLAS/MKL threads
to one before numeric imports. Single-core behavior is algorithmic/thread-env
controlled; CPU affinity is not pinned on macOS.

Workers compare actual site permutations, independent Burnside orbit counts,
and (for retained records) orbit-size coverage. Canonical representatives are
normalized outside timing to the same minimum convention, sorted with full
stabilizer/orbit-size metadata, checked unique, and SHA-256 hashed. CASM's raw
maximum convention and ID order can differ from ours; normalized orbit records
must agree. Verification cost is excluded on both sides.

Output JSONL is checkpointed after each attempt. The JSON and CSV contain raw
records and median summaries. Existing checkpoints are **resumed**, not rerun.
For a new measurement, use a fresh path:

```sh
.venv/bin/python experiment/casm_scale_compare.py --output validation/new_scale_run.json
# A smaller pilot:
.venv/bin/python experiment/casm_scale_compare.py --case fjh 1 3 --reps 1 --output validation/new_pilot.json
```

The watchdog defaults to 180 s total process duration and 1536 MiB sampled RSS.
These include outside-timer setup/verification. Time or memory truncations are
recorded as censored outcomes, not exhaustive times or zero counts. Matched
speedups are emitted only when all attempted runs for both sides are exhaustive.

The committed 2026-11-16 study completed all 180 attempts exhaustively. See
`validation/CASM_SCALE_REPORT_2026-11-16.md` for results, version details, source
fingerprints and interpretation. Existing paper/old comparison CSV are preserved.
