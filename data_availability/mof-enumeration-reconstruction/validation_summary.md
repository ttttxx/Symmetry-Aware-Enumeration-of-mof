# Reconstruction Audit Summary

Audit date: 2026-10-10

Protocol: `tables/reconstruction_audit_protocol_seeded_2026-10-10.json`

## Scope and sampling

- `ffc`: all 42 symmetry-distinct assignments in the configured building-unit pool.
- `fjh`: 30 of 394 orbit representatives, sampled with a fixed seed at 10 assignments per High, Med, and Low stabilizer stratum.
- Total: 72 reconstructed structures. The protocol file records the exact selected IDs and seed used for each topology.

## Results

| Topology | Audited | Assembled | Sanitized CIF | Topology recovered | Clash-free | Closure pairs | Closure RMSD (Å) | Pairs ≤0.8 Å |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `ffc` | 42 | 42/42 | 42/42 | 42/42 | 42/42 | 1,008 | 0.515 | 83.3% |
| `fjh` | 30 | 30/30 | 30/30 | 30/30 | 30/30 | 2,880 | 0.332 | 100.0% |
| Combined | 72 | 72/72 | 72/72 | 72/72 | 72/72 | 3,888 | 0.388 | 95.7% |

Every CIF passed the sanitizer's clean status, reloaded successfully, recovered its target quotient-graph topology, and received no short-contact flag. The closure statistics summarize node–linker connection-site pairs recorded during reconstruction. DFT is not part of this structural audit.

The earlier 27-structure audit summary is retained in the CSV files without the `seeded_2026-10-10` suffix. The manuscript's current results use the fixed-seed audit summarized above.
