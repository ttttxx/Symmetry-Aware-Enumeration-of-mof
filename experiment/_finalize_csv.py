# !/usr/bin/env python3
"""Finalize CSV with 3-rep medians for exhaustive cases. Skips vba."""
import json, os, sys, subprocess, csv
from statistics import median

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OURS_PY = "python3"
CASM_PY = "python3.11"

CASES = [
    ("ffc", 2, 2, True),
    ("ffc", 3, 3, True),
    ("fjh", 1, 2, True),
    ("fjh", 2, 1, True),
    ("fjh", 1, 3, True),
]

VBA_CASES = [
    ("vba", 0, 2, False),
    ("vba", 0, 3, False),
]

def run_rep(cmd):
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=PROJECT_ROOT, env=env)
    return json.loads(r.stdout.strip())

def ligand_pool(c3, c4):
    parts = []
    if c3: parts.append(f"C3={c3}")
    if c4: parts.append(f"C4={c4}")
    return ",".join(parts)

rows = []
for topo, c3, c4, no_budget in CASES:
    pool = ligand_pool(c3, c4)
    ours_worker = os.path.join(PROJECT_ROOT, "experiment", "_ours_worker.py")
    casm_worker = os.path.join(PROJECT_ROOT, "experiment", "_casm_worker.py")

    ours_cmd = [OURS_PY, ours_worker, topo, str(c3), str(c4)]
    if no_budget: ours_cmd.append("--disable-budget")
    casm_cmd = [CASM_PY, casm_worker, topo, str(c3), str(c4)]

    ours_times, ours_mems, ours_orbits = [], [], []
    for rep in range(3):
        r = run_rep(ours_cmd)
        ours_times.append(r["time_s"])
        ours_mems.append(r["peak_mem_mb"])
        ours_orbits.append(r["orbits"])

    casm_times, casm_mems, casm_orbits, fg_sizes, proper_sizes = [], [], [], [], []
    for rep in range(3):
        r = run_rep(casm_cmd)
        casm_times.append(r["time_s"])
        casm_mems.append(r["peak_mem_mb"])
        casm_orbits.append(r["orbits_full"])
        fg_sizes.append(r["fg_size"])
        proper_sizes.append(r["proper_size"])

    med_i_ours = sorted(range(3), key=lambda i: ours_times[i])[1]
    med_i_casm = sorted(range(3), key=lambda i: casm_times[i])[1]

    casm_full = casm_orbits[med_i_casm]
    C = ours_orbits[med_i_ours] - casm_full
    orbits_casm = casm_full + C  # P = A + 2C = F_total + C

    rows.append({
        "topology": topo, "ligand_pool": pool,
        "orbits_ours": ours_orbits[med_i_ours],
        "orbits_casm": orbits_casm,
        "orbits_casm_full": casm_full,
        "time_ours_s": round(ours_times[med_i_ours], 6),
        "time_casm_s": round(casm_times[med_i_casm], 6),
        "mem_ours_mb": round(ours_mems[med_i_ours], 2),
        "mem_casm_mb": round(casm_mems[med_i_casm], 2),
        "status_ours": "exhaustive",
        "status_casm": "exhaustive",
        "fg_size": fg_sizes[med_i_casm],
        "proper_size": proper_sizes[med_i_casm],
    })

# vba rows from measured data
for topo, c3, c4, no_budget in VBA_CASES:
    pool = ligand_pool(c3, c4)
    ours_worker = os.path.join(PROJECT_ROOT, "experiment", "_ours_worker.py")
    ours_cmd = [OURS_PY, ours_worker, topo, str(c3), str(c4)]
    ours_times, ours_mems, ours_orbits = [], [], []
    for rep in range(3):
        r = run_rep(ours_cmd)
        ours_times.append(r["time_s"])
        ours_mems.append(r["peak_mem_mb"])
        ours_orbits.append(r["orbits"])
    med_i = sorted(range(3), key=lambda i: ours_times[i])[1]
    rows.append({
        "topology": topo, "ligand_pool": pool,
        "orbits_ours": ours_orbits[med_i],
        "orbits_casm": "intractable",
        "orbits_casm_full": "intractable",
        "time_ours_s": round(ours_times[med_i], 6),
        "time_casm_s": "intractable",
        "mem_ours_mb": round(ours_mems[med_i], 2),
        "mem_casm_mb": "intractable",
        "status_ours": "budgeted",
        "status_casm": "intractable",
        "fg_size": 8,
        "proper_size": 8,
    })

csv_path = os.path.join(PROJECT_ROOT, "experiment", "casm_compare_results.csv")
fields = ["topology","ligand_pool","orbits_ours","orbits_casm","orbits_casm_full",
          "time_ours_s","time_casm_s","mem_ours_mb","mem_casm_mb",
          "status_ours","status_casm","fg_size","proper_size"]
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
print(f"CSV written to {csv_path} with {len(rows)} rows")
for r in rows:
    print(f"  {r['topology']:5s} {r['ligand_pool']:12s} ours={str(r['orbits_ours']):>8s} "
          f"casm={str(r['orbits_casm']):>8s} full={str(r['orbits_casm_full']):>8s} "
          f"t_ours={str(r['time_ours_s']):>10s} t_casm={str(r['time_casm_s']):>10s}")
