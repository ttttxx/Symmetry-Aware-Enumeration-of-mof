from __future__ import annotations

import json
import os
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from co2_n2_screening.raspa import parse_gcmc_labels


PROJECT = Path("/Users/jing/PyCharmMiscProject/mof")
CO2_N2 = Path("/Users/jing/PyCharmMiscProject/CO2_N2")
RUN_ROOT = PROJECT / "results/co2n2_repeat_run_2026-10-19/replicate_002"
PROTOCOL = PROJECT / "results/co2n2_matched_ordering_expansion_2026-10-10/protocol.json"
SOURCE_JOBS = (
    CO2_N2
    / "outputs/co2n2_ffc_screening/batch_002_clean_pool_5000/gcmc_jobs"
)
SIMULATE = CO2_N2 / ".conda-tools/bin/simulate"
RASPA_DIR = CO2_N2 / ".conda-tools"


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    groups = protocol["matched_groups"]
    candidate_ids = sorted(
        {cid for group in groups.values() for cid in group["candidate_ids"]}
    )
    RUN_ROOT.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, object]] = []
    jobs: list[dict[str, str]] = []

    for cid in candidate_ids:
        source = SOURCE_JOBS / cid
        target = RUN_ROOT / cid
        target.mkdir()
        for phase in ("adsorption", "desorption"):
            src_phase = source / phase
            dst_phase = target / phase
            dst_phase.mkdir()
            for filename in ("simulation.input", f"{cid}.cif"):
                source_file = src_phase / filename
                target_file = dst_phase / filename
                if not source_file.is_file():
                    raise FileNotFoundError(source_file)
                target_file.write_bytes(source_file.read_bytes())
        jobs.append({"candidate_id": cid, "job_dir": str(target)})

    manifest = {
        "replicate": "independent_replicate_002",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_source": str(PROTOCOL),
        "source_jobs": str(SOURCE_JOBS),
        "simulation_binary": str(SIMULATE),
        "raspa_dir": str(RASPA_DIR),
        "candidate_count": len(candidate_ids),
        "independent_phase_runs": 2 * len(candidate_ids),
        "input_seed_override": None,
        "seed_policy": "RASPA default; record the emitted seed from each phase output",
        "candidate_ids": candidate_ids,
    }
    (RUN_ROOT / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    env = os.environ.copy()
    env["RASPA_DIR"] = str(RASPA_DIR)
    for job in jobs:
        cid = job["candidate_id"]
        start_candidate = time.monotonic()
        phase_records = []
        for phase in ("adsorption", "desorption"):
            phase_dir = RUN_ROOT / cid / phase
            phase_start = time.monotonic()
            result = subprocess.run(
                [str(SIMULATE), "simulation.input"],
                cwd=phase_dir,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            (phase_dir / "raspa_stdout.log").write_text(
                result.stdout, encoding="utf-8"
            )
            (phase_dir / "raspa_stderr.log").write_text(
                result.stderr, encoding="utf-8"
            )
            data_files = sorted(
                (phase_dir / "Output/System_0").glob("*.data")
            )
            output_text = (
                data_files[-1].read_text(encoding="utf-8", errors="replace")
                if data_files
                else ""
            )
            seed_match = re.search(r"Random number seed:\s*(\d+)", output_text)
            phase_record = {
                "candidate_id": cid,
                "phase": phase,
                "returncode": result.returncode,
                "elapsed_s": round(time.monotonic() - phase_start, 3),
                "output_file": str(data_files[-1]) if data_files else "",
                "random_seed": int(seed_match.group(1)) if seed_match else None,
                "finished": "Simulation finished" in output_text,
            }
            phase_records.append(phase_record)
            records.append(phase_record)
            (RUN_ROOT / "phase_runs.json").write_text(
                json.dumps(records, indent=2), encoding="utf-8"
            )
        status = "succeeded" if all(
            row["returncode"] == 0 and row["finished"] for row in phase_records
        ) else "failed"
        print(
            f"{cid}: {status}; elapsed={time.monotonic() - start_candidate:.1f}s",
            flush=True,
        )
        jobs_record = {**job, "gcmc_status": status, "gcmc_error": ""}
        job_path = RUN_ROOT / "jobs.csv"
        pd.DataFrame(jobs + []).to_csv(
            job_path, index=False
        ) if not job_path.exists() else None
        if status == "failed":
            (RUN_ROOT / "failed_candidates.txt").open("a").write(cid + "\n")

    job_table = pd.DataFrame(
        [
            {
                "candidate_id": job["candidate_id"],
                "job_dir": job["job_dir"],
                "gcmc_status": "succeeded"
                if all(
                    row["returncode"] == 0 and row["finished"]
                    for row in records
                    if row["candidate_id"] == job["candidate_id"]
                )
                else "failed",
                "gcmc_error": "",
            }
            for job in jobs
        ]
    )
    job_table.to_csv(RUN_ROOT / "jobs.csv", index=False)
    labels = parse_gcmc_labels(job_table, RUN_ROOT / "replicate_002_labels.csv")
    labels.to_csv(RUN_ROOT / "replicate_002_labels.csv", index=False)
    print(f"Wrote {len(labels)} candidate labels to {RUN_ROOT}", flush=True)


if __name__ == "__main__":
    main()
