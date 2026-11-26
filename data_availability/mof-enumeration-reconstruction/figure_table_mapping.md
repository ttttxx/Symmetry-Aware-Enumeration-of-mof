# Figure and Table Mapping

## Main Manuscript

| Manuscript item | Package files |
|---|---|
| Fig. 1 workflow schematic | Generated manuscript artwork; no tabular source data in this package |
| Fig. 2 reconstruction gallery, `ffc` panel | `cifs/representative/fig2_ffc_representative.cif` |
| Fig. 2 reconstruction gallery, `fjh` panel | `cifs/representative/fig2_fjh_representative.cif` |
| Fig. 2 reconstruction gallery, `vba` panel | `cifs/representative/fig2_vba_representative.cif` |
| Reconstruction audit result in the Results text | `tables/reconstruction_records_seeded_2026-10-10.csv`, `tables/reconstruction_by_topology_seeded_2026-10-10.csv`, `tables/reconstruction_overall_seeded_2026-10-10.csv`, `tables/reconstruction_audit_protocol_seeded_2026-10-10.json`, `cifs/reconstruction_audit_seeded_2026-10-10/` |
| Matched CO2/N2 ordering check table (candidate-level means of five independent GCMC runs) | `tables/ffc_co2n2_ordering_check.csv`, `tables/ffc_co2n2_gcmc_run1.csv`–`tables/ffc_co2n2_gcmc_run5.csv`, `tables/co2n2_five_run_mean_protocol.json`, `tables/co2n2_five_run_seed_manifest.csv` |
| Bounded `vba` exact-$\rho$ output table | `reproducibility/unified_20261123/MATCH_UNIFIED_VBA_2026-11-23.json` (per-layer medians, counts, vectors and ordered digests) |
| X-site fusion audit in Supplementary Information | `tables/x_site_fusion_*.csv`, `tables/x_site_fusion_protocol.md`, and the corresponding audited CIFs in `cifs/reconstruction_audit_ffc_fjh/` |
| Enumeration / CASM comparison table | `reproducibility/unified_20261123/MATCH_UNIFIED_REFRESH_2026-11-23.json` and `.csv` (current ABI-4, matched Gamma and records) |
| Single-point DFT screening figures | `tables/single_point_dft_source_data.csv`, `cifs/dft_screening_2026-10-20/` |

## Notes

The representative, reconstruction-audit, and DFT input CIF files are pre-relaxation starting models selected for figure inspection or the stated analyses. The current quantitative reconstruction audit includes every `ffc` orbit representative and a fixed-seed, symmetry-stratified `fjh` sample. The older 27-structure files remain in the package as provenance. The CO2/N2 ordering-check CSV is a compact extract from an existing screening run and is used for the matched fixed-composition screening example.

## Current supplementary tables

| Manuscript item | Package files |
|---|---|
| Matched-group records, Gamma/F counts, time ranges/RSS, native core and CASM count-only tables | `reproducibility/unified_20261123/MATCH_UNIFIED_REFRESH_2026-11-23.json` |
| G1–G9 node compositions | `tables/ffc_co2n2_group_compositions.csv`, checked against all 30 `ffc_co2n2_ordering_check.csv` records |
| Enumeration implementation and exact input version | `reproducibility/unified_20261123/MATCH_UNIFIED_MANIFEST_2026-11-23.json`, `snapshot/` |

## Approved manuscript wording (2026-11-24)

The current manuscript source is `reproducibility/manuscript_20261124/main_MATCH.tex`. All numerical enumeration tables retain the 2026-11-23 measurement sources above; the dated measured snapshot remains unchanged.
