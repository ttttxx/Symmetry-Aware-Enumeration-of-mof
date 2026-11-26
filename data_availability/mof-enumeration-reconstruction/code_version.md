# Code Version Notes

Prepared date: 2026-08-28

Working source repository observed during package preparation:

- Remote: `https://github.com/ttttxx/Symmetry-Aware-Enumeration-of-mof.git`
- Source commit: `e35ee3451c19d76d7b659b84ecb11ee116a2d21c`
- Manuscript and data-availability edits were prepared in the local working tree after this commit.

The manuscript Data availability statement points to the intended public data/code repository:

`https://github.com/Buct454/mof-enumeration-reconstruction.git`

Before submission, copy or push the package contents into that public repository and make sure the public commit hash or release tag is recorded here.

On 2026-10-20, the local package was extended with the CASM comparison table, 53 converged single-point DFT values and their input CIFs, the nine added Zeo++ descriptors and protocol, and the X-site fusion cutoff/validation records. The package manifest now lists individual files, and the figure/table map includes these datasets. Upload and public-link verification remain pending; no remote repository files were changed in this preparation step.

## Reproducible reconstruction audit update

On 2026-10-10, the structural audit was rerun with `experiment/geometric_reconstruction_summary.py`. It covered all 42 `ffc` orbit representatives and 30 of 394 `fjh` representatives, selected with a fixed seed and stratified by stabilizer class. The command, seed derivation, class counts, selected IDs, and per-structure outputs are recorded in `tables/reconstruction_audit_protocol_seeded_2026-10-10.json`. This audit runs assembly and structural checks; it does not invoke DFT.

## Bounded `vba` exact-rho output experiment

On 2026-10-10, `experiment/vba_rho_strata_benchmark.py` compared the legacy DFS prefix with 20,000-result ceilings at each of the five stabilizer orders for the `vba` group (`|Gamma|=16`). Both the two-unit Cu/B pool and the three-unit Cu/B/Cd pool used the existing 100,000 total output ceiling and a 60 s fixed-point-task timeout. Legacy prefixes were additionally rerun at the exact total returned by the quota policy, enabling matched-output layer comparisons. Enumeration logs, protocol, summaries, and compressed representative lists are in `results/vba_rho_stratification_20261010/`; the manuscript table data are in `tables/vba_rho_stratification_20261010.csv`. This experiment measures output allocation and enumeration cost, not structural quality or screening performance. The prior traversal-order dataset is retained as `tables/vba_order_sensitivity_20260914.csv`.

## Current ABI-4 enumeration replacement (2026-11-23)

All manuscript enumeration timings now come from the current-source rerun in `reproducibility/unified_20261123/`. Exact source/library/input hashes are in its manifest, with the measured files preserved in `snapshot/`. Both matched-group and bounded-output experiments use the same lexical template-pool convention. The default CLI/API backend is now auto (C++ preferred, optimized Python fallback); benchmarks explicitly required C++. The new local files have not been pushed or verified as public by this task.

## Approved wording update (2026-11-24)

Fourteen approved review items were applied to the manuscript, with the source in `reproducibility/manuscript_20261124/`. Mathematical equations, measured code, input fingerprints and numerical tables are unchanged. Linker scope now explicitly uses the same template on every edge, matching the experiments. The Cu role entry omits unreported open-shell probes. No remote upload or publication was performed.
