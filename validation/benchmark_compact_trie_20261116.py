"""Fresh-process native hash/compact-Trie comparison, same DFS and outputs."""
import argparse
import contextlib
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import statistics
import struct
import subprocess
import sys
import time

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CASES = ["ffc_3_3", "ffc_5_8", "fjh_1_3", "fjh_1_4", "lil_4_3", "swap18_partial"]


def worker(case, store):
    import sympy.combinatorics
    from enumerator_native import HybridMOFEnumerator
    from experiment.casm_scale_common import pools_for, orbit_count
    from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
    from symmetry import SymmetryEngine
    from validation.benchmark_prefix_trie_20261116 import InterruptedSwapEnumerator
    with contextlib.redirect_stdout(sys.stderr):
        if case == "swap18_partial":
            n = 18
            perms = [tuple(range(n)), (1, 0) + tuple(range(2, n))]
            pools, degrees, options = {0: ["A", "B"]}, [0] * n, {}
            cls = InterruptedSwapEnumerator
        else:
            topology, c3, c4 = case.split("_")
            graph = enrich_graph_data_with_topology(load_graph(str(ROOT/"data"/topology/f"{topology}_graph.json")))
            pools = pools_for(int(c3), int(c4))
            degrees = [node["degree"] for node in graph["nodes"]]
            perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
            options = dict(user_selected_ligands=pools, node_degrees=degrees,
                           topology_requirements=calculate_topology_requirements(graph))
            n, cls = len(degrees), HybridMOFEnumerator
        labels = sorted({label for pool in pools.values() for label in pool})
        expected = orbit_count(perms, degrees, pools)
        naive = math.prod(len(pools[d]) for d in degrees)
        instance = cls(perms, n, backend="cpp", native_store=store, prefix_trie=True,
                       verbose=False, max_total_results=None, identity_orbit_limit=float("inf"), **options)
        started = time.perf_counter()
        records = instance.enumerate(labels)
        elapsed = time.perf_counter() - started
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak_mib = peak / (1048576 if sys.platform == "darwin" else 1024)
        assert instance.active_backend == "cpp" and instance.active_store == "cpp_" + store
        assert len(records) == expected
        assert sum(record["orbit_size"] for record in records) == naive
        # Stream the exact sorted output (including original insertion IDs),
        # avoiding an additional giant JSON string after peak RSS measurement.
        digest = hashlib.sha256()
        codes = {label:i for i,label in enumerate(labels)}
        for record in records:
            digest.update(struct.pack("<QII", record["id"], record["stabilizer_size"], record["orbit_size"]))
            digest.update(bytes(codes[label] for label in record["coloring"]))
        sample = records[::max(1,len(records)//100)]
        assert all(instance.classes[instance._trie_lookup(tuple(record["coloring"]))] is record for record in sample)
        return dict(case=case,store=store,sites=n,group_size=len(perms),naive=naive,
                    orbits=len(records),expected_burnside=expected,time_s=elapsed,
                    peak_rss_mib=peak_mib,output_sha256=digest.hexdigest(),
                    kernel_s=sum(t["kernel_s"] for t in instance.native_task_reports),
                    task_stats=instance.task_dedup_stats,prefix_stats=instance.prefix_trie_stats,
                    index_stats=instance.result_index_stats,
                    python_hash_entries=len(instance._canonical_index))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", choices=CASES)
    parser.add_argument("--store", choices=("hash", "trie"))
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker, args.store)))
        return
    sources = ["native/enumeration_core.cpp", "enumerator_native.py", "enumerator.py",
               "enumerator_task_dedup.py", "enumerator_casm_inspired.py", "enumerator_storage.py",
               "search_prefix_trie.py", "experiment/casm_scale_common.py", "main.py",
               "validation/benchmark_prefix_trie_20261116.py", "validation/benchmark_compact_trie_20261116.py",
               "native/libmof_enumeration" + (".dylib" if sys.platform == "darwin" else ".so")]
    hashes = {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources}
    records = []
    for case in CASES:
        for repetition in range(args.repetitions):
            for store in (("hash", "trie") if repetition % 2 == 0 else ("trie", "hash")):
                result = subprocess.run([sys.executable, __file__, "--worker", case, "--store", store],
                                        capture_output=True, text=True, timeout=180)
                if result.returncode:
                    raise RuntimeError(result.stderr)
                row = json.loads(result.stdout)
                row["rep"] = repetition
                records.append(row)
                print(case, repetition, store, round(row["time_s"],4), round(row["peak_rss_mib"],1),flush=True)
    summary = []
    for case in CASES:
        selected = [row for row in records if row["case"] == case]
        assert len({row["output_sha256"] for row in selected}) == 1
        assert len({row["task_stats"]["search_nodes"] for row in selected}) == 1
        assert len({row["task_stats"]["leaf_visits"] for row in selected}) == 1
        result = dict(case=case,orbits=selected[0]["orbits"],all_records_and_ids_match=True)
        for store in ("hash", "trie"):
            rows = [row for row in selected if row["store"] == store]
            result[store] = dict(median_s=statistics.median(r["time_s"] for r in rows),
                                 kernel_median_s=statistics.median(r["kernel_s"] for r in rows),
                                 peak_rss_median_mib=statistics.median(r["peak_rss_mib"] for r in rows),
                                 times_s=[r["time_s"] for r in rows])
        summary.append(result)
    assert all(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest for name,digest in hashes.items())
    output = dict(python=sys.version,platform=platform.platform(),native_abi=3,repetitions=args.repetitions,
                  timing_scope="fresh process enumerate; preparation/imports/validation outside; RSS before digest",
                  prefix_trie=True,native_kernel="dfs",source_sha256=hashes,summary=summary,records=records)
    path=ROOT/"validation/COMPACT_TRIE_BENCHMARK_2026-11-16.json"
    path.write_text(json.dumps(output,indent=2)+"\n")
    print(path)


if __name__ == "__main__":
    main()
