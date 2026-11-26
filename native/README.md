# Optional C++ search core

Python retains graph/group analysis, subgroup scheduling, fixed-point-domain
exclusion compilation, result dictionaries, ordering and CIF construction.
C++17 performs a task's DFS, canonical-prefix comparison, stabilizer counting
and exact complete-key checks using a compact array Trie by default. Results cross the ctypes boundary
in batches, once per subgroup task, rather than once per search node.

Build on macOS/Linux with an installed C++17 compiler (no extra binding package):

```sh
.venv/bin/python tools/build_enumeration_core.py
```

The generated shared library is platform-specific and git-ignored. Rebuild on
the target machine. Importing the Python wrapper never runs a compiler.

Use the existing CLI pipeline (its existing topology/ligand settings still apply):

```sh
.venv/bin/python main.py --enumeration-backend cpp
.venv/bin/python main.py --enumeration-backend auto
.venv/bin/python main.py --enumeration-backend python
```

`cpp` requires the native library; `auto` falls back to the optimized Python
byte-key DFS if loading fails or inputs exceed native limits. `python` forces
that Python backend. Omitting the option selects `auto`: C++ is preferred, with optimized Python fallback. Select `legacy` explicitly for the original implementation.
The `run_enumeration` Python API accepts the same `enumeration_backend` option.
Benchmark JSON records requested and active backends. Exact-rho quotas now share the optimized search core; see the ABI-4 quota section below.

The initial native implementation supports 1–512 sites and at most 256 encoded
labels. Input label order, canonical convention, result fields, sorting, total
output cap and identity-task new-result budget match the Python implementation.
Coverage masks use multiple 64-bit words, so the number of covered regions is
not limited to 64. The comparison hash mode compares full keys on collisions.

The unified algorithm, state definitions, correctness argument and pseudocode
are in [docs/subgroup_dfs_trie_algorithm.md](../docs/subgroup_dfs_trie_algorithm.md).
The normal C++ result index is a path-compressed byte Trie: a contiguous node
array, 16 bytes per node, with integer child/sibling indices. Edge labels refer
to slices of existing result keys. The implementation allocates neither a map
per node nor a second copy of each edge label. C++ uses no full-key hash table
in Trie mode, and Python retains the exported flat index instead of a full-key
hash mirror. Sorted output IDs differ from the stable insertion indices used
inside the Trie; post-run lookup returns the latter, matching the existing API.

```sh
# Full algorithm: compact representative index and completed-prefix reuse
.venv/bin/python main.py --enumeration-backend cpp --native-store trie --prefix-trie
# Same native DFS with the former hash index, for comparison
.venv/bin/python main.py --enumeration-backend cpp --native-store hash --prefix-trie
```

The constructor accepts `native_store="trie"` (default) or `"hash"`.
`run_enumeration` accepts the same optional argument. This setting controls C++
storage; Python fallback still uses its existing exact byte/tuple hash index.
Benchmark JSON records the active store and index statistics. `node_bytes`
counts logical Trie nodes; `capacity_bytes` includes spare vector capacity.
Neither includes common result-key/output storage. The current native ABI is
**4**: rebuild the platform library. Trie node/record indices are 32-bit, checked
before overflow. See `validation/COMPACT_TRIE_REVIEW_2026-11-16.md` for measurements.

Task deadlines bound native search using a steady clock; argument preparation
is deducted before entry. Like the Python path, this is a cooperative cutoff.
Bulk result transfer/materialization can finish after the search deadline.
Timed-out tasks are conservatively excluded from whole-domain coverage
certificates; fixed-time output sets can differ because execution speed changes.
The shared library releases the GIL through ctypes during the search call.

Tests and same-language DFS/counter comparison:

```sh
.venv/bin/python -m pytest -q tests/test_native_enumerator.py tests/test_enumerator_storage.py tests/test_counter_enumerator.py tests/test_subgroup_task_dedup.py tests/test_identity_budget.py
.venv/bin/python validation/benchmark_hybrid_20261113.py
```

The `native_kernel="counter"` Python constructor option is for the comparison
experiment. Normal CLI use keeps subgroup DFS.

## Completed-prefix Trie (opt-in)

```sh
.venv/bin/python tools/build_enumeration_core.py
.venv/bin/python main.py --enumeration-backend cpp --prefix-trie
```

The Python backend supports the same flag. The API accepts `prefix_trie=True`
and `prefix_trie_max_nodes=50000`. Prefix caching remains opt-in; this is a
separate search-coverage index from the complete-representative Trie. This
prefix feature originally used ABI 2; the current complete-key array
implementation originally required ABI 3; exact-order filtering now requires ABI 4.

A closed Trie prefix certifies that its source task's canonical DFS subtree
has been exhausted. On interruption, the active DFS frontier is compressed into
closed sibling prefixes; the unfinished path remains open. A later task can
skip a subtree only when every completion satisfies the source partition's
equal-color constraints, its allowed label domains fit inside the source
domains, and a closed prefix matches. Prefix existence alone never permits a
prune. Queries are scheduled at depths where the relevant values and equalities
first become decidable. No per-node Python callbacks are added to C++ search.

Fully completed tasks continue to use the existing whole-domain exclusions and
allocate no prefix nodes. The additional opportunity is completed portions of
interrupted tasks. A global node cap bounds cached Trie nodes (not total
process memory); reaching it drops coverage certificates and leaves search
correct. The cache resets on each enumeration. C++ frees its cache with the
native context. `prefix_trie_stats` reports queries, hits, `early_hits` (before
leaf depth), closed prefixes, node counts and budget skips. CLI benchmark JSON
includes these statistics. `search_nodes` counts entered DFS frames; it excludes
assignments rejected by canonical comparison before recursive entry.

The feature stays opt-in: complete exhaustive runs already have whole-task
coverage, and small saved regions can cost more to query than they save. Exact
Exact-rho quotas scope coverage certificates to one stabilizer layer; the comparison counter does not use completed-prefix queries.
Time budgets remain cooperative; creating the bounded frontier snapshot can
finish after a task's search deadline.

```sh
.venv/bin/python -m pytest -q tests/test_search_prefix_trie.py
.venv/bin/python validation/benchmark_prefix_trie_20261116.py
```

See `validation/PREFIX_TRIE_REVIEW_2026-11-16.md` for the on/off results and
limitations. The earlier CASM measurements used ABI 1; their measured sources
and local binary were archived under `validation/snapshots/casm_scale_20261116`
before editing. Those times do not measure this prefix-Trie change.


## Exact-stabilizer quotas (ABI 4)

```sh
.venv/bin/python tools/build_enumeration_core.py
.venv/bin/python main.py --enumeration-backend cpp --native-store trie --prefix-trie --rho-quota-per-level 20000
```

`HybridRhoQuotaMOFEnumerator` shares the original quota schedule through
`RhoQuotaScheduleMixin`. Tasks are ordered by decreasing subgroup order and
lexicographic partition order within each layer; labels retain supplied order.
The fused Python/native search accepts a full leaf only when its stabilizer
order equals the current target. Rejected orders never enter the result index
or consume a quota. The task's added-result budget is the remaining layer quota,
and the global output cap is enforced independently.

Completed-task equalities and completed-prefix certificates are scoped to the
current exact stabilizer layer. Changing layers clears both kinds of coverage
certificate while preserving the global complete-key index and all outputs.
A timed-out task supplies only completed-prefix certificates, never whole-task
coverage. Reaching a quota is reported as `quota_reached`, not exhaustion.
`python` forces optimized integer DFS; `auto` can fall back to it; `cpp` requires
the rebuilt library. The default CLI/API backend is `auto`; `legacy` remains available explicitly.

The original quota scheduler's positional deadline call was also corrected to
pass `deadline=` explicitly; previously it occupied the identity-budget slot.
Historical results/CSV are retained. No previous timing is relabelled as ABI 4.

Tests: `tests/test_native_rho_quota.py`. Real-pool verification and fresh-process
measurements: `validation/benchmark_rho_native_20261123.py` and
`validation/RHO_NATIVE_VBA_2026-11-23.json`.
