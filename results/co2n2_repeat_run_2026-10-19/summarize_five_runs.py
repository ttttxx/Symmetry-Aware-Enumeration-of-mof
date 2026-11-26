from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/jing/PyCharmMiscProject/mof")
RUN_ROOT = ROOT / "results/co2n2_repeat_run_2026-10-19"
PACKAGE = ROOT / "data_availability/mof-enumeration-reconstruction"
TABLES = PACKAGE / "tables"
PROTOCOL = ROOT / "results/co2n2_matched_ordering_expansion_2026-10-10/protocol.json"
OLD_TABLE = TABLES / "ffc_co2n2_ordering_check.csv"
DIRECT_METRICS = [
    "co2_uptake_mmolg",
    "n2_uptake_mmolg",
    "co2_desorption_uptake_mmolg",
    "co2_working_capacity_mmolg",
    "co2_n2_selectivity",
    "regenerability",
]


def main() -> None:
    labels: dict[int, pd.DataFrame] = {}
    for rep_no in range(1, 6):
        rep_dir = RUN_ROOT / f"replicate_{rep_no:03d}"
        if rep_no == 1:
            src = RUN_ROOT / "replicate_002/replicate_001_reparsed_labels.csv"
        elif rep_no == 2:
            src = RUN_ROOT / "replicate_002/replicate_002_labels.csv"
        else:
            src = rep_dir / f"replicate_{rep_no:03d}_labels.csv"
        frame = pd.read_csv(src)
        if len(frame) != 30 or not frame["gcmc_status"].eq("succeeded").all():
            raise ValueError(f"Incomplete or failed labels in {src}")
        frame = frame.set_index("candidate_id")
        if rep_no == 1:
            expected_ids = set(frame.index)
        elif set(frame.index) != expected_ids:
            raise ValueError(f"Candidate set mismatch in replicate {rep_no}")
        labels[rep_no] = frame
        package_copy = TABLES / f"ffc_co2n2_gcmc_run{rep_no}.csv"
        shutil.copy2(src, package_copy)

    metrics = pd.concat(
        {rep: frame[DIRECT_METRICS] for rep, frame in labels.items()},
        names=["replicate", "candidate_id"],
    ).reset_index()
    metrics["aps"] = metrics["co2_n2_selectivity"] * metrics["co2_working_capacity_mmolg"]
    summary = metrics.groupby("candidate_id").agg(
        **{f"{metric}_mean": (metric, "mean") for metric in DIRECT_METRICS},
        **{f"{metric}_sd": (metric, "std") for metric in DIRECT_METRICS},
        aps_run_mean=("aps", "mean"),
        aps_sd=("aps", "std"),
    )
    summary["aps_mean"] = summary["co2_n2_selectivity_mean"] * summary["co2_working_capacity_mmolg_mean"]
    summary["aps_cv_percent"] = summary["aps_sd"] / summary["aps_mean"] * 100
    summary = summary.reset_index()
    summary["gcmc_runs"] = 5

    base = pd.read_csv(OLD_TABLE).drop(
        columns=["source_run1_labels", "source_run2_labels"], errors="ignore"
    )
    base = base.drop(columns=[c for c in base.columns if c.endswith("_sd") or c.endswith("_mean") or c == "aps_cv_percent"], errors="ignore")
    base = base.drop(columns=["gcmc_runs"], errors="ignore")
    combined = base.merge(summary, on="candidate_id", validate="one_to_one")
    for metric in DIRECT_METRICS:
        combined[metric] = combined[f"{metric}_mean"]
    combined["aps"] = combined["aps_mean"]
    combined["gcmc_runs"] = 5
    for rep_no in range(1, 6):
        combined[f"source_run{rep_no}_labels"] = f"tables/ffc_co2n2_gcmc_run{rep_no}.csv"
    combined.to_csv(OLD_TABLE, index=False, float_format="%.10g")
    combined.to_csv(RUN_ROOT / "candidate_metrics_five_run_summary.csv", index=False, float_format="%.10g")

    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    group_rows = []
    winners: dict[str, list[str]] = {}
    for group, info in protocol["matched_groups"].items():
        block = combined[combined["candidate_id"].isin(info["candidate_ids"])].copy()
        row = {"group": group, "rho": info["rho"], "n": len(block)}
        for metric in ["pld_angstrom", "lcd_angstrom", "co2_working_capacity_mmolg", "co2_n2_selectivity", "regenerability", "aps"]:
            row[f"{metric}_min"] = block[metric].min()
            row[f"{metric}_max"] = block[metric].max()
        row["aps_mean_sd"] = block["aps_sd"].mean()
        row["aps_max_sd"] = block["aps_sd"].max()
        group_rows.append(row)
        winners[group] = []
        for rep_no, frame in labels.items():
            ids = info["candidate_ids"]
            run = frame.loc[ids].copy()
            run["aps"] = run["co2_n2_selectivity"] * run["co2_working_capacity_mmolg"]
            winners[group].append(str(run["aps"].idxmax()))
    group_summary = pd.DataFrame(group_rows)
    group_summary.to_csv(RUN_ROOT / "matched_group_ranges_five_run_mean.csv", index=False, float_format="%.10g")
    winner_summary = {
        group: {
            "winning_candidate_each_run": values,
            "modal_winner": max(set(values), key=values.count),
            "modal_winner_count_of_5": values.count(max(set(values), key=values.count)),
        }
        for group, values in winners.items()
    }
    (RUN_ROOT / "top_aps_winner_stability.json").write_text(json.dumps(winner_summary, indent=2) + "\n")

    phase_counts = {}
    new_seeds = []
    all_seeds = []
    for rep_no in range(2, 6):
        records = json.loads((RUN_ROOT / f"replicate_{rep_no:03d}/phase_runs.json").read_text())
        phase_counts[str(rep_no)] = len(records)
        all_seeds.extend(int(row["random_seed"]) for row in records)
        if rep_no >= 3:
            new_seeds.extend(int(row["random_seed"]) for row in records)
    original_seed = []
    seed_rows = []
    for candidate in expected_ids:
        for phase in ("adsorption", "desorption"):
            for path in (Path("/Users/jing/PyCharmMiscProject/CO2_N2/outputs/co2n2_ffc_screening/batch_002_clean_pool_5000/gcmc_jobs") / candidate / phase / "Output/System_0").glob("*.data"):
                import re
                match = re.search(r"Random number seed:\s*(\d+)", path.read_text(errors="replace"))
                if match:
                    seed = int(match.group(1))
                    original_seed.append(seed)
                    seed_rows.append({"replicate": 1, "candidate_id": candidate, "phase": phase, "random_seed": seed, "status": "succeeded"})
    for rep_no in range(2, 6):
        records = json.loads((RUN_ROOT / f"replicate_{rep_no:03d}/phase_runs.json").read_text())
        for row in records:
            seed_rows.append({
                "replicate": rep_no,
                "candidate_id": row["candidate_id"],
                "phase": row["phase"],
                "random_seed": int(row["random_seed"]),
                "status": "succeeded" if row["returncode"] == 0 and row["finished"] else "failed",
            })
    seed_manifest = pd.DataFrame(seed_rows).sort_values(["replicate", "candidate_id", "phase"])
    if len(seed_manifest) != 300 or seed_manifest["random_seed"].duplicated().any():
        raise ValueError("Expected 300 successful phase runs with globally unique seeds")
    seed_manifest.to_csv(TABLES / "co2n2_five_run_seed_manifest.csv", index=False)
    previous_seeds = set(original_seed)
    previous_seeds.update(int(row["random_seed"]) for row in json.loads((RUN_ROOT / "replicate_002/phase_runs.json").read_text()))
    all_seeds = original_seed + all_seeds
    protocol_summary = {
        "independent_runs_per_candidate": 5,
        "candidate_count": 30,
        "phase_runs_per_candidate_per_replicate": 2,
        "new_replicates": [3, 4, 5],
        "phase_run_counts_by_replicate": phase_counts,
        "additional_replicate_phase_run_total": sum(len(json.loads((RUN_ROOT / f"replicate_{rep_no:03d}/phase_runs.json").read_text())) for rep_no in (3, 4, 5)),
        "new_seed_count": len(new_seeds),
        "first_run_seed_count": len(original_seed),
        "new_seeds_unique": len(new_seeds) == len(set(new_seeds)),
        "all_five_run_seeds_unique": len(all_seeds) == len(set(all_seeds)),
        "new_seeds_disjoint_from_first_two_runs": not previous_seeds.intersection(new_seeds),
        "direct_metric_mean_rule": "Each direct GCMC metric is the arithmetic mean across five independent runs.",
        "aps_rule": "APS is calculated from the five-run mean selectivity and five-run mean CO2 working capacity; run-to-run SD is calculated from per-run APS values.",
        "uptake_unit": "mmol/g from RASPA Average loading absolute [mol/kg framework]",
        "cycles_per_phase": 10000,
        "initialization_cycles_per_phase": 5000,
    }
    (TABLES / "co2n2_five_run_mean_protocol.json").write_text(json.dumps(protocol_summary, indent=2) + "\n", encoding="utf-8")
    print("Per-candidate APS CV (%):")
    print(summary["aps_cv_percent"].describe(percentiles=[0.5, 0.9]).to_string())
    print("\nGroup APS ranges and average run SD:")
    print(group_summary[["group", "aps_min", "aps_max", "aps_mean_sd", "aps_max_sd"]].to_string(index=False))
    print("\nTop APS winner stability:")
    print(json.dumps({k: v["modal_winner_count_of_5"] for k, v in winner_summary.items()}, indent=2))
    print("\nProtocol checks:")
    print(json.dumps(protocol_summary, indent=2))


if __name__ == "__main__":
    main()
