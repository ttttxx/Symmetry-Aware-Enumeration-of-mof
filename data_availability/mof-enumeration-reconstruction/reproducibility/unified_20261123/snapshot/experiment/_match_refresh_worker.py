"""Single-process, single-repetition worker for matched-group scale benchmarks."""
import argparse
import collections
import contextlib
import gc
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import sys
import time

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import sympy.combinatorics
import libcasm.configuration as config
import libcasm.enumerate as enum
from libcasm.enumerate._enumerate import ConfigEnumAllOccupationsBase
from enumerator_native import HybridMOFEnumerator
from enumerator import MOFEnumerator
from experiment.casm_scale_common import prepare, record_digest


def run(topology,c3,c4,side,progress):
    def status(phase,**extra):
        if progress:
            Path(progress).write_text(json.dumps(dict(phase=phase,**extra)))
    status("preparing")
    inputs = prepare(topology,c3,c4)
    meta = inputs["meta"]
    degrees,pools = inputs["degrees"],inputs["pools"]
    labels = sorted({label for pool in pools.values() for label in pool})
    records = None
    instance = None
    if side in ("ours", "legacy"):
        cls = HybridMOFEnumerator if side == "ours" else MOFEnumerator
        extra = dict(backend="cpp", native_store="trie", prefix_trie=True) if side == "ours" else {}
        instance = cls(inputs["perms"],len(degrees),verbose=False,**extra,
            topology_requirements=inputs["requirements"],user_selected_ligands=pools,
            node_degrees=degrees,max_total_results=None,identity_orbit_limit=float("inf"),task_time_limit=None)
    gc.collect()
    status("measuring")
    t0 = time.perf_counter()
    if side in ("ours", "legacy"):
        records = instance.enumerate(labels)
        if side == "ours":
            assert instance.active_backend == "cpp" and instance.active_store == "cpp_trie"
        count = len(records)
    elif side == "casm_full":
        count = sum(1 for _ in enum.ConfigEnumAllOccupations(inputs["prim"]).by_supercell_list([inputs["cell"]]))
    else:
        base = ConfigEnumAllOccupationsBase(config.Configuration(inputs["cell"]),set(range(len(degrees))))
        proper = inputs["proper"]
        count,visited = 0,0
        by_site = [pools[degree] for degree in degrees]
        if side == "casm_records":
            records = []
        while base.is_valid():
            current = base.value()
            visited += 1
            if config.is_canonical_configuration(current,subgroup=proper):
                if records is not None:
                    stabilizer = len(config.make_invariant_subgroup(current,group=proper))
                    coloring = [pool[value] for pool,value in zip(by_site,current.occupation.tolist())]
                    records.append(dict(id=count,coloring=coloring,stabilizer_size=stabilizer,
                                         orbit_size=len(proper)//stabilizer))
                count += 1
            base.advance()
        assert visited == meta["naive"]
        if records is not None:
            records.sort(key=lambda row:row["stabilizer_size"],reverse=True)
            for index,row in enumerate(records):
                row["id"] = index
    elapsed = time.perf_counter()-t0
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1048576 if sys.platform=="darwin" else 1024)
    status("verifying",elapsed_s=elapsed,orbits=count)
    expected = meta["expected_full_orbits"] if side == "casm_full" else meta["expected_orbits"]
    assert count == expected,(side,count,expected)
    result = dict(**meta,side=side,time_s=elapsed,peak_rss_mib=peak,orbits=count,
                  status="exhaustive",python=sys.version.split()[0],
                  casm_version=importlib.metadata.version("libcasm-configuration"))
    if records is not None:
        assert sum(row["orbit_size"] for row in records) == meta["naive"]
        result["record_digest"] = record_digest(records,labels,inputs["perms"])
        result["stabilizer_histogram"] = dict(collections.Counter(row["stabilizer_size"] for row in records))
    if side == "ours":
        result["prefix_stats"] = instance.prefix_trie_stats
        result["index_stats"] = instance.result_index_stats
        result["backend"] = instance.active_backend
        result["store"] = instance.active_store
        result["native_kernel_s"] = sum(row["kernel_s"] for row in instance.native_task_reports)
        result["task_stats"] = instance.task_dedup_stats
    status("complete",elapsed_s=elapsed,orbits=count)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("topology")
    parser.add_argument("c3",type=int)
    parser.add_argument("c4",type=int)
    parser.add_argument("side",choices=["ours","legacy","casm_records","casm_count","casm_full"])
    parser.add_argument("--progress")
    args = parser.parse_args()
    with contextlib.redirect_stdout(sys.stderr):
        result = run(args.topology,args.c3,args.c4,args.side,args.progress)
    print(json.dumps(result))
