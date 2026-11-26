from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/jing/PyCharmMiscProject/mof")
BASE = ROOT / "results/co2n2_repeat_run_2026-10-19/replicate_002"
RUN_ROOT = ROOT / "results/co2n2_repeat_run_2026-10-19"
SIMULATE = Path("/Users/jing/PyCharmMiscProject/CO2_N2/.conda-tools/bin/simulate")
RASPA_DIR = SIMULATE.parent.parent
ORIGINAL_JOBS = Path("/Users/jing/PyCharmMiscProject/CO2_N2/outputs/co2n2_ffc_screening/batch_002_clean_pool_5000/gcmc_jobs")
LOADING = re.compile(
    r"Average loading absolute\s+\[mol/kg framework\]\s+([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)"
)
SEED = re.compile(r"Random number seed:\s*(\d+)")


def parse_phase(output: Path) -> tuple[dict[str, float], int]:
    text = output.read_text(encoding="utf-8", errors="replace")
    component = ""
    loads: dict[str, float] = {}
    for line in text.splitlines():
        m = re.search(r"Component\s+\d+\s+[\[(]?([A-Za-z0-9_+-]+)[\])] ?", line)
        if m:
            component = m.group(1).upper()
        elif "CO2" in line.upper():
            component = "CO2"
        elif "N2" in line.upper():
            component = "N2"
        lm = LOADING.search(line)
        if lm and component:
            loads[component] = float(lm.group(1))
    seed_match = SEED.search(text)
    if seed_match is None:
        raise RuntimeError(f"No RASPA random seed found in {output}")
    if not {"CO2", "N2"}.issubset(loads):
        raise RuntimeError(f"Missing component loading in {output}: {loads}")
    return loads, int(seed_match.group(1))


def one_phase(candidate: str, phase: str, replicate_dir: Path) -> tuple[dict[str, float], dict[str, object]]:
    source = BASE / candidate / phase
    target = replicate_dir / candidate / phase
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source / "simulation.input", target / "simulation.input")
    cifs = list(source.glob("*.cif"))
    if len(cifs) != 1:
        raise RuntimeError(f"Expected one CIF in {source}, found {len(cifs)}")
    shutil.copy2(cifs[0], target / cifs[0].name)
    env = os.environ.copy()
    env["RASPA_DIR"] = str(RASPA_DIR)
    started = time.monotonic()
    proc = subprocess.run([str(SIMULATE)], cwd=target, env=env, text=True, capture_output=True)
    (target / "raspa_stdout.log").write_text(proc.stdout, encoding="utf-8")
    (target / "raspa_stderr.log").write_text(proc.stderr, encoding="utf-8")
    outputs = list((target / "Output/System_0").glob("*.data"))
    if proc.returncode != 0 or len(outputs) != 1:
        raise RuntimeError(f"RASPA failed for {candidate}/{phase}: rc={proc.returncode}, outputs={outputs}")
    values, seed = parse_phase(outputs[0])
    return values, {
        "candidate_id": candidate,
        "phase": phase,
        "returncode": proc.returncode,
        "elapsed_s": round(time.monotonic() - started, 3),
        "output_file": str(outputs[0]),
        "random_seed": seed,
        "finished": True,
    }


def main() -> None:
    if not SIMULATE.is_file():
        raise FileNotFoundError(SIMULATE)
    candidates = pd.read_csv(BASE / "replicate_002_labels.csv")["candidate_id"].tolist()
    for rep_no in (3, 4, 5):
        rep_name = f"replicate_{rep_no:03d}"
        rep_dir = RUN_ROOT / rep_name
        rep_dir.mkdir(parents=True, exist_ok=True)
        started_utc = datetime.now(timezone.utc).isoformat()
        records: list[dict[str, object]] = []
        labels: list[dict[str, object]] = []
        seeds: list[int] = []
        print(f"Starting {rep_name}: {len(candidates)} candidates, 60 phase runs", flush=True)
        for i, candidate in enumerate(candidates, start=1):
            phase_values = {}
            for phase in ("adsorption", "desorption"):
                values, record = one_phase(candidate, phase, rep_dir)
                records.append(record)
                phase_values[phase] = values
                seeds.append(int(record["random_seed"]))
            co2_ads = phase_values["adsorption"]["CO2"]
            n2_ads = phase_values["adsorption"]["N2"]
            co2_des = phase_values["desorption"]["CO2"]
            wc = co2_ads - co2_des
            selectivity = (co2_ads / n2_ads) / (0.15 / 0.85) if n2_ads else float("nan")
            regenerability = wc / co2_ads if co2_ads else float("nan")
            labels.append({
                "candidate_id": candidate,
                "co2_uptake_mmolg": co2_ads,
                "n2_uptake_mmolg": n2_ads,
                "co2_desorption_uptake_mmolg": co2_des,
                "co2_working_capacity_mmolg": wc,
                "co2_n2_selectivity": selectivity,
                "regenerability": regenerability,
                "aps": selectivity * wc,
                "gcmc_status": "succeeded",
                "gcmc_error": "",
            })
            if i % 5 == 0:
                print(f"{rep_name}: {i}/{len(candidates)} candidates complete", flush=True)
        if len(set(seeds)) != len(seeds):
            raise RuntimeError(f"Duplicate seeds in {rep_name}")
        prior_seeds = {
            int(x["random_seed"]) for x in json.loads((BASE / "phase_runs.json").read_text())
        }
        for candidate in candidates:
            for phase in ("adsorption", "desorption"):
                for output in (ORIGINAL_JOBS / candidate / phase / "Output/System_0").glob("*.data"):
                    match = SEED.search(output.read_text(encoding="utf-8", errors="replace"))
                    if match:
                        prior_seeds.add(int(match.group(1)))
        if prior_seeds.intersection(seeds):
            raise RuntimeError(f"Seeds overlap an earlier run in {rep_name}")
        pd.DataFrame(labels).to_csv(rep_dir / f"{rep_name}_labels.csv", index=False)
        (rep_dir / "phase_runs.json").write_text(json.dumps(records, indent=2) + "\n")
        manifest = {
            "replicate": f"independent_{rep_name}",
            "created_utc": started_utc,
            "candidate_count": len(candidates),
            "independent_phase_runs": len(records),
            "seed_policy": "RASPA default; emitted seed recorded; seeds checked unique and disjoint from replicate_002",
            "simulation_binary": str(SIMULATE),
        }
        (rep_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        print(f"Completed {rep_name}: 60/60 phase runs; all seeds unique", flush=True)


if __name__ == "__main__":
    main()
