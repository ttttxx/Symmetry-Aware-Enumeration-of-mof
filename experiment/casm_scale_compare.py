"""Fresh-process matched-group gradient experiment, checkpointed after every run."""
import argparse
import csv
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from experiment.casm_scale_common import CASES
import psutil

SIDES = ["ours","casm_records","casm_count","casm_full"]


def launch(python,case,side,timeout_s,memory_mib):
    with tempfile.TemporaryDirectory(prefix="mof-casm-scale-") as directory:
        progress = Path(directory)/"progress.json"
        with (Path(directory)/"stdout").open("w+") as out,(Path(directory)/"stderr").open("w+") as err:
            command = [python,str(ROOT/"experiment"/"_casm_scale_worker.py"),*map(str,case),side,
                       "--progress",str(progress)]
            process = subprocess.Popen(command,stdout=out,stderr=err,cwd=ROOT)
            monitor = psutil.Process(process.pid)
            started = time.monotonic()
            observed_peak = 0
            reason = None
            while process.poll() is None:
                try:
                    rss = monitor.memory_info().rss/1048576
                    observed_peak = max(observed_peak,rss)
                except psutil.Error:
                    pass
                if observed_peak > memory_mib:
                    reason = "memory_limit"
                elif time.monotonic()-started > timeout_s:
                    reason = "wall_timeout"
                if reason:
                    process.kill()
                    process.wait()
                    break
                time.sleep(.1)
            out.seek(0); err.seek(0)
            stdout,stderr = out.read(),err.read()
            if reason or process.returncode:
                phase = json.loads(progress.read_text()) if progress.exists() else {}
                return dict(topology=case[0],c3=case[1],c4=case[2],side=side,
                    status=reason or "error",observed_peak_rss_mib=observed_peak,
                    progress=phase,stderr=stderr[-2000:],returncode=process.returncode)
            result = json.loads(stdout)
            result["observed_peak_rss_mib"] = observed_peak
            return result


def summarize(records,output):
    groups = {}
    for record in records:
        key = (record["topology"],record["c3"],record["c4"])
        groups.setdefault(key,{}).setdefault(record["side"],[]).append(record)
    rows = []
    for case,sides in groups.items():
        row = dict(topology=case[0],c3=case[1],c4=case[2])
        digests = set()
        digest_sides = set()
        for side,runs in sides.items():
            complete = [r for r in runs if r["status"] == "exhaustive"]
            row[side+"_complete_reps"] = len(complete)
            row[side+"_status"] = "exhaustive" if len(complete)==len(runs) else ";".join(sorted({r["status"] for r in runs}))
            if complete:
                row["naive"] = complete[0]["naive"]
                row["orbits"] = complete[0]["expected_orbits"]
                row["group_size"] = complete[0]["group_size"]
                row[side+"_median_s"] = statistics.median(r["time_s"] for r in complete)
                row[side+"_median_rss_mib"] = statistics.median(r["peak_rss_mib"] for r in complete)
                digests.update(r["record_digest"] for r in complete if "record_digest" in r)
                if any("record_digest" in r for r in complete):
                    digest_sides.add(side)
        if len(digests)>1:
            raise AssertionError(f"Normalized complete records differ: {case}")
        row["matched_record_digests"] = len(digests)==1 and {"ours","casm_records"} <= digest_sides
        if all(row.get(side+"_status")=="exhaustive" for side in ["ours","casm_records"]):
            row["matched_records_speedup"] = row["casm_records_median_s"]/row["ours_median_s"]
        if all(row.get(side+"_status")=="exhaustive" for side in ["ours","casm_count"]):
            row["vs_casm_count_speedup"] = row["casm_count_median_s"]/row["ours_median_s"]
        rows.append(row)
    output.write_text(json.dumps(dict(records=records,summary=rows),indent=2))
    if rows:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with output.with_suffix(".csv").open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    return rows


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python",default="python3.11")
    parser.add_argument("--reps",type=int,default=3)
    parser.add_argument("--output",type=Path,default=ROOT/"validation"/"CASM_SCALE_2026-11-16.json")
    parser.add_argument("--timeout",type=float,default=180)
    parser.add_argument("--memory-mib",type=float,default=1536)
    parser.add_argument("--case",nargs=3,action="append",metavar=("TOPOLOGY","C3","C4"))
    args=parser.parse_args()
    cases=[(t,int(c3),int(c4)) for t,c3,c4 in args.case] if args.case else CASES
    args.output.parent.mkdir(parents=True,exist_ok=True)
    checkpoint=args.output.with_suffix(".jsonl")
    records=[json.loads(line) for line in checkpoint.read_text().splitlines()] if checkpoint.exists() else []
    completed={(r["topology"],r["c3"],r["c4"],r["side"],r["rep"]) for r in records}
    for case_index,case in enumerate(cases):
        for rep in range(args.reps):
            offset=(case_index+rep)%len(SIDES)
            for side in SIDES[offset:]+SIDES[:offset]:
                if (*case,side,rep) in completed:
                    continue
                result=launch(args.python,case,side,args.timeout,args.memory_mib)
                result.update(rep=rep,watchdog_timeout_s=args.timeout,watchdog_memory_mib=args.memory_mib)
                records.append(result)
                with checkpoint.open("a") as stream:
                    stream.write(json.dumps(result)+"\n")
                rows=summarize(records,args.output)
                elapsed=result.get("time_s")
                print(case,side,f"rep {rep+1}",result["status"],
                      f"{elapsed:.4f}s" if elapsed is not None else result.get("progress"),flush=True)
                if result["status"]=="error":
                    raise RuntimeError(result.get("stderr"))
        current=next(row for row in rows if (row["topology"],row["c3"],row["c4"])==case)
        print("CASE",json.dumps(current),flush=True)
