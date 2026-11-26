"""Validate the completed scale study and generate a compact result table."""
from pathlib import Path
import hashlib
import json
import statistics
import csv
import argparse

ROOT=Path(__file__).resolve().parents[1]


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=ROOT,
                        help="Measured source tree, e.g. the archived pre-prefix snapshot")
    parser.add_argument("--check-only", action="store_true", help="Validate without rewriting results")
    args=parser.parse_args()
    manifest=json.loads((ROOT/"validation/CASM_SCALE_MANIFEST_2026-11-16.json").read_text())
    data=json.loads((ROOT/"validation/CASM_SCALE_2026-11-16.json").read_text())
    expected={(*case,side,rep) for case in manifest["cases"] for side in manifest["modes"]
              for rep in range(manifest["repetitions"])}
    observed={(r["topology"],r["c3"],r["c4"],r["side"],r["rep"]) for r in data["records"]}
    assert observed==expected,("Experiment incomplete",len(observed),len(expected))
    assert len(data["records"])==len(expected),"Repeated checkpoint rows"
    assert all(hashlib.sha256((args.source_root/name).read_bytes()).hexdigest()==digest
               for name,digest in manifest["source_sha256"].items()),"Measured source or native binary changed"
    rows=[]
    for row in data["summary"]:
        case=(row["topology"],row["c3"],row["c4"])
        selected=[r for r in data["records"] if (r["topology"],r["c3"],r["c4"])==case]
        complete=all(r["status"]=="exhaustive" for r in selected)
        result=dict(row,all_runs_exhaustive=complete)
        if complete:
            assert row["matched_record_digests"]
            assert all(r["orbits"]==r["expected_full_orbits" if r["side"]=="casm_full" else "expected_orbits"] for r in selected)
            assert len({r["record_digest"] for r in selected if "record_digest" in r})==1
            assert len({json.dumps(r["stabilizer_histogram"],sort_keys=True) for r in selected if "stabilizer_histogram" in r})==1
            ours=[r for r in selected if r["side"]=="ours"]
            result["casm_full_group_orbits"]=selected[0]["expected_full_orbits"]
            result["casm_full_group_size"]=selected[0]["casm_full_size"]
            result["ours_native_kernel_median_s"]=statistics.median(r["native_kernel_s"] for r in ours)
            result["ours_non_kernel_remainder_s"]=row["ours_median_s"]-result["ours_native_kernel_median_s"]
            for side in manifest["modes"]:
                times=[r["time_s"] for r in selected if r["side"]==side]
                result[side+"_range_s"]=[min(times),max(times)]
        rows.append(result)
    output=dict(cases=len(rows),attempts=len(data["records"]),
                exhaustive_attempts=sum(r["status"]=="exhaustive" for r in data["records"]),
                source_hashes_unchanged=True,rows=rows)
    if args.check_only:
        print(json.dumps(dict(source_root=str(args.source_root.resolve()),
                              **{k:v for k,v in output.items() if k!="rows"})))
        raise SystemExit(0)
    (ROOT/"validation/CASM_SCALE_ANALYSIS_2026-11-16.json").write_text(json.dumps(output,indent=2))
    with (ROOT/"validation/CASM_SCALE_ANALYSIS_2026-11-16.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows({key:json.dumps(value) if isinstance(value,(list,dict)) else value
                          for key,value in row.items()} for row in rows)
    lines=["| 拓扑 | C3/C4 | 原始组合 | Γ轨道 | 本文/s | CASM同结果/s | 加速比 | 对CASM只计数加速比 |",
           "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        if row["all_runs_exhaustive"]:
            lines.append(f'| {row["topology"]} | {row["c3"]}/{row["c4"]} | {row["naive"]:,} | {row["orbits"]:,} | '
                         f'{row["ours_median_s"]:.6f} | {row["casm_records_median_s"]:.6f} | '
                         f'{row["matched_records_speedup"]:.2f}× | {row["vs_casm_count_speedup"]:.2f}× |')
        else:
            lines.append(f'| {row["topology"]} | {row["c3"]}/{row["c4"]} | — | — | 资源截断 | — | — | — |')
    table="\n".join(lines)+"\n"
    (ROOT/"validation/CASM_SCALE_TABLE_2026-11-16.md").write_text(table)
    print(json.dumps({k:v for k,v in output.items() if k!="rows"}))
    print(table)
