# CASM scale study measured sources

Copied and SHA-256 verified against `CASM_SCALE_MANIFEST_2026-11-16.json`
before implementing the completed-prefix Trie. Contains the 13 measured
source/library paths, including the local arm64 ABI-1 shared library. This is
a provenance snapshot, not a standalone runnable distribution. Do not update
these copies to match later source changes.

Validate the historical 180-run result set against these measured files:

```sh
.venv/bin/python validation/analyze_casm_scale_20261116.py \
  --source-root validation/snapshots/casm_scale_20261116 --check-only
```

The live wrapper/core now use ABI 2. The CASM study's recorded times remain
measurements of the archived version. Its manifest and raw results are unchanged.
