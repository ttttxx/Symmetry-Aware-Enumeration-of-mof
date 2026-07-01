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

Recorded a server-side DFT interruption pattern for future resume/debugging:

- A traceback ending with `KeyboardInterrupt: received signal 15; stopped CP2K` comes from `dft_validation/cp2k_runner.py`, which converts external `SIGTERM` into a controlled stop of the whole CP2K process group.
- Interpretation: the nohup Python driver did not simply keep running in the background after that message; the current batch was interrupted and needs explicit status checking or `--resume`.
- `tail -f "$LOG"` failing with `cannot open ''` was traced to an empty `LOG` environment variable in the current shell, most likely because the original variable assignment was made in an earlier shell session.
- Practical follow-up on the server: check `ps -fp "$(cat run_parallel.pid)"` if the pid file exists, otherwise `ps aux | grep '[r]un_dft_validation.py'`; if nothing is alive, restart with `--resume` and record the new log filename explicitly.

Reviewed the current `ffc` evidence against likely reviewer criticism and tightened the manuscript accordingly:

- Confirmed from the latest `selected_structures.json` that `ffc` currently has 21 converged single-point calculations across 8 composition groups.
- Composition-wise minima fall in `Med/High` symmetry for all 8 groups (`High=6`, `Med=2`, `Low=0`), which supports a limited candidate-retention claim.
- Also confirmed that this dataset is too small and too imbalanced across compositions to defend as a strong statistical symmetry--energy trend.
- Updated `paper/main.tex` to reflect this more carefully: `ffc` is now framed as a limited qualitative check, not a universal trend claim; the unrelaxed single-point protocol is described more explicitly as a screening-stage comparison on a common assembled-geometry manifold.
- Recorded a design warning for future experiments: removing composition-controlled minima markers from the `ffc` pooled figure, or replacing them with only global minima, would weaken the argument rather than strengthen it.
- Recorded a code warning: adding `3c_Cu_imadazoline_Ch` to `ffc.atom_choices` would change the chemistry pool, but by itself does not guarantee a stronger validation dataset and may worsen composition confounding.

Assessed `fjh` runtime expectations and acceleration options:

- A wall time of roughly 3 hours per `fjh` single-point on the current CPU-only shared server is plausible for the present setup: ~760-atom cells, diagonalization SCF, and PBE+D3 with `500/60 Ry`.
- If the scientific settings must remain unchanged, the safest way to reduce total calendar time is throughput parallelism rather than trying to make each individual `fjh` job dramatically faster.
- The first practical knob to test is lower OpenMP per job with more concurrent `fjh` jobs, e.g. `CP2K_OMP=16` and `DFT_MAX_WORKERS["fjh"]=2`, provided memory stays within budget and the server is not crowded.
- More aggressive accelerations such as lower cutoffs, smaller basis, or weaker SCF settings would change the protocol and should not be adopted for the main paper results without an explicit methods decision.
