# Worklog

## 2026-06-30

Created persistent Codex project memory under `codex/`.

Captured current MOF enumeration project state:

- Core graph/group-theory enumeration model and Symmetry Index definition.
- Existing main modules and their responsibilities.
- Current DFT validation strategy using composition-controlled comparisons.
- Rationale for moving from small `ffc` validation to richer `fjh` topology.
- CP2K single-point DFT method, server paths, and runbook.
- Paper wording direction: qualitative validation and pre-screening heuristic, not a stability guarantee.

Next time, start by reading `codex/README.md`, `codex/PROJECT_STATE.md`, and `codex/NEXT_ACTIONS.md`.

Updated DFT validation execution and paper figure wiring:

- Added batch-level CP2K concurrency in `dft_validation/run_dft_validation.py`.
- Kept CP2K itself as `.ssmp` OpenMP; concurrency now means multiple independent single-point jobs in parallel.
- Added resource knobs in `dft_validation/config.py`: `DFT_PARALLEL`, `DFT_MAX_WORKERS`, and `DFT_JOB_MEM_GB`.
- Kept conservative defaults for the shared Lab207 server: `ffc=4` concurrent jobs, `fjh=2` concurrent jobs, `CP2K_OMP=16`.
- Clarified that analysis outputs separate topology figures: `rho_energy_ffc.png` and `rho_energy_fjh.png`.
- Updated `paper/main.tex` to use separate `ffc` and `fjh` panels instead of the old `physical_validation_relative.png` placeholder.
- Fixed a stale config comment so production DFT consistently says `500 / 60 Ry`.

Operational caution: before restarting with the new parallel runner on the server, stop the old nohup Python/CP2K process tree to avoid duplicate runs for the same structures.

Refined validation strategy after questioning whether `fjh` 2/rho is sufficient:

- Increased recommended `fjh` sampling to 2 composition groups and up to 4 structures per distinct `rho`.
- Added `--harvest-only` to parse existing `cp2k_run.out` files into `selected_structures.json` without launching CP2K.
- Added `--summary` to print a compact energy table for checking whether parsed DFT energies are physically plausible.
- Interpretation: `fjh` should be the main evidence; `ffc` should not be discarded, but treated as small-topology sanity check and used immediately through harvesting.

Updated paper validation framing and figure:

- Reframed the claim from monotonic `rho`-energy decrease to low-energy candidates being retained/enriched in medium- and high-symmetry regions.
- Generated `paper/ffc_pooled_symmetry_energy.png` from the downloaded `dft_vs_symmetry.csv`.
- Used fixed symmetry bins: Low `0 < rho <= 0.25`, Medium `0.25 < rho <= 0.75`, High `0.75 < rho <= 1`.
- Inserted the pooled `ffc` figure into `paper/main.tex`.
- Added English paragraphs plus `%：中文：` comments explaining the non-relaxed single-point protocol as an enumeration-time pre-screening score, not relaxed thermodynamic stability.

Refined figure strategy again:

- Removed the duplicate composition-controlled `fig:symmetry_stability` from the main DFT validation narrative.
- Added paired pooled figures for `ffc` and `fjh`: `paper/ffc_pooled_symmetry_energy.png` and `paper/fjh_pooled_symmetry_energy.png`.
- Rewrote the validation text in a more defensible form: the claim is candidate retention in medium/high symmetry, not monotonic `rho`-energy decrease.
- Added a caption caveat that the current `fjh` panel is preliminary until all selected `fjh` symmetry strata finish.
- Recompiled `paper/main.tex` successfully with `latexmk`.

Tightened the paper defense around pooled DFT figures:

- Confirmed `fig:symmetry_stability` is no longer present in `paper/main.tex`; the active DFT figure label is `fig:pooled_symmetry_energy`.
- Clarified that the paired `ffc`/`fjh` pooled figures are descriptive topology-level evidence.
- Made the primary visual criterion explicit: composition-wise lowest-energy open markers, not a global pooled slope.
- Strengthened the single-point DFT caveat: energies are geometry-consistent screening scores, not relaxed formation energies or a full stability ranking.
- Recompiled `paper/main.tex` successfully after the wording update.
