# X-site fusion audit data

The audit uses the 27 fixed-geometry structures listed in the manuscript: 12 `ffc` and 15 `fjh` reconstructions sampled across symmetry strata with base seed `20260703`. Their CIFs are included under `cifs/reconstruction_audit_ffc_fjh/`.

- `x_site_fusion_model_sensitivity.csv` records each model's cutoff-dependent matched pairs and topology-linked versus non-topological merges over 0.4--1.2 Å.
- `x_site_fusion_model_sensitivity_summary.csv` gives the corresponding per-cutoff aggregate.
- `x_site_fusion_threshold_scan.csv` records the broader fixed-geometry scan over 0.4--2.0 Å, including cutoff-dependent missed expected pairs and false pairs.
- `x_site_fusion_topology_linked_pairs.csv` lists template-linked node--linker pairs and their bond distances.
- `x_site_fusion_full_validation_summary.csv` reports complete reconstruction and CIF-reload checks at 0.8 and 1.2 Å.

The broader scan evaluates geometric matching. Full reconstruction and CIF validation were performed at 0.8 and 1.2 Å, as summarized in the manuscript. These records support the fixed-geometry cutoff analysis; they do not represent relaxed structures.
