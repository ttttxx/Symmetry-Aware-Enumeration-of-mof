# Zeo++ descriptor correction, 2026-11-10

The `-res` numeric columns are Di (largest included sphere diameter), Df (largest free sphere diameter), and Dif (largest included sphere diameter along the free-sphere path). The corrected mapping is LCD = Di and PLD = Df.

The original screening adapter assigned Di to PLD and Df to LCD. All 30 matched candidates were checked against their original `.res` files. The current `ffc_co2n2_ordering_check.csv` and nine-row `targeted_zeopp_descriptors.csv` have been re-exported with the corrected mapping. The manuscript pore ranges and current five-run group summaries use the corrected fields. The raw geometry calculations were retained; no new Zeo++ or GCMC simulation was required.

GCMC loadings, selectivity, working capacity, working-capacity fraction, APS, run-to-run statistics, candidate membership, and composition groups were retained. Across the nine groups, within-group spans are 0.03–1.02 Å for PLD and 0.18–0.92 Å for LCD.

Audit evidence and the re-export script are maintained in the MOF project under `codex/validation/match_pld_lcd_raw_audit_20261110.csv` and `codex/validation/correct_match_pld_lcd_20261110.py`. Dated earlier intermediate snapshots retain their original values and should not be used as corrected descriptor tables. The upstream CO2_N2 adapter has also been corrected; its older batch exports require explicit re-export before reuse.

Zeo++ output format: https://github.com/mharanczyk/zeoplusplus/blob/master/documentation/README
