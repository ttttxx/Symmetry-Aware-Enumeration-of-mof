# Next Actions

Last updated: 2026-06-30

## Immediate

- Check current DFT run status on `huangbo@Lab207` before changing or deleting any DFT output directories.
- Harvest already completed CP2K outputs before deciding to rerun anything:
  `python -u dft_validation/run_dft_validation.py --resume --harvest-only --summary`
- Before deploying the parallel runner, stop the old serial nohup Python job and confirm no old CP2K processes remain; otherwise the same structures may be run twice.
- If DFT jobs are running, monitor `run_all.log` and the active `cp2k_run.out`.
- When new DFT results complete, run or inspect `dft_validation/analyze.py` outputs grouped by topology and composition.
- Until `fjh` finishes, keep the paper claim on `ffc` explicitly limited to a qualitative retention check, not a statistical symmetry--energy trend.

## Analysis

- For each fixed composition, report Spearman rank correlation between `rho` and energy per atom.
- Check whether the lowest-energy structure in each composition group falls in the maximum-`rho` subset; use continuous `rho` as the primary symmetry variable.
- Treat mixed-composition trends cautiously; do not pool all structures as evidence for symmetry-stability correlation.
- Use separate figures for topology-level evidence: `rho_energy_ffc.png` and `rho_energy_fjh.png`.
- Treat `fjh` as the stronger validation evidence because it has more nodes and five distinct `rho` levels within fixed composition; use `ffc` mainly as a smaller sanity-check topology.
- Prefer expanding `fjh` to 2 composition groups and up to 4 structures per `rho` before making the paper claim too strong.
- Keep the paper claim as candidate-retention/enrichment in medium/high symmetry, not monotonic `rho`-energy decrease.
- Keep the `ffc` open-circle markers unless a replacement figure still preserves composition-controlled minima explicitly; removing them and highlighting only global minima would weaken the logic.
- In the main paper, keep `ffc_pooled_symmetry_energy.png` and `fjh_pooled_symmetry_energy.png` as parallel panels under `fig:pooled_symmetry_energy`; do not reintroduce the old duplicate `fig:symmetry_stability` main-text figure unless it is moved to supplementary material with a distinct purpose.
- After more `fjh` DFT jobs finish, regenerate `paper/fjh_pooled_symmetry_energy.png` from the updated CSV/JSON and remove or soften the "preliminary" caveat if all selected strata are complete.
- If reviewers may question non-relaxed single points, consider a small relaxed sanity check for representative `fjh` points after the main single-point screening is complete.
- Do not assume that adding `3c_Cu_imadazoline_Ch` to `ffc.atom_choices` will automatically strengthen the paper. It changes the chemistry pool and may worsen composition confounding unless the validation design is revised intentionally.

## Code Consistency

- Keep `dft_validation/select_structures.py` aligned with the composition-controlled stratified sampling design.
- Keep `dft_validation/run_dft_validation.py` topology handling generic; avoid reintroducing hard-coded `ffc` assumptions.
- Keep `dft_validation/cp2k_runner.py` aligned with the current `.ssmp` OpenMP CP2K setup and single-point method.
- Keep `dft_validation/run_dft_validation.py` batch concurrency conservative on the shared 128 GB server, especially for `fjh`.

## Paper Consistency

- Keep `paper/main.tex` consistent with actual DFT settings:
  - PBE+D3.
  - `DZVP-MOLOPT-PBE-GTH`.
  - `500 / 60 Ry`.
  - Single-point `ENERGY`.
  - Composition-controlled stratified sampling.
- Avoid claiming that higher symmetry guarantees lower energy.
- Present symmetry-prioritized enumeration as an a priori pre-screening heuristic.
