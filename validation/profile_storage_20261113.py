"""Fresh-process storage benchmark plus separately instrumented phase/memory runs.

Default: five cases x three indexes x three independent timing processes, then
largest-case phase and tracemalloc runs. Each imports SymPy before the timer.
No tracing/instrumentation is active in normal timing/RSS measurements.
"""
import builtins
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
BACKENDS = ["trie", "tuple", "bytes"]
CASES = [("ffc", 2, 2), ("ffc", 3, 3), ("fjh", 1, 2), ("fjh", 2, 1), ("fjh", 1, 3)]


def worker(backend, topology, c3, c4, mode):
    import gc
    import math
    import sympy.combinatorics
    from experiment._ours_worker import build_ligand_dict
    from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
    from symmetry import SymmetryEngine
    import enumerator
    from enumerator_casm_inspired import CasmInspiredMOFEnumerator
    from enumerator_storage import HashStorageMOFEnumerator

    graph = enrich_graph_data_with_topology(load_graph(str(ROOT / "data" / topology / f"{topology}_graph.json")))
    perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
    pools = build_ligand_dict(c3, c4)
    labels = sorted({label for pool in pools.values() for label in pool})
    kwargs = dict(verbose=False, topology_requirements=calculate_topology_requirements(graph),
                  user_selected_ligands=pools, node_degrees=[node["degree"] for node in graph["nodes"]],
                  max_total_results=None, identity_orbit_limit=float("inf"))
    base = CasmInspiredMOFEnumerator if backend == "trie" else HashStorageMOFEnumerator
    if backend != "trie":
        kwargs["storage_key"] = backend
    phases = {}
    cls = base
    if mode == "phase":
        cls = type("PhaseMeasured", (base,), {})
        for name in ["enumerate_subgroups", "_compile_exclusions", "_search_compiled_task",
                     "_store_canonical_leaf", "_trie_lookup", "_trie_insert", "_verify_burnside"]:
            original = getattr(base, name)

            def make_wrapper(name, original):
                def wrapped(self, *args, **kw):
                    start = time.perf_counter()
                    result = original(self, *args, **kw)
                    phases[name] = phases.get(name, 0.) + time.perf_counter()-start
                    return result
                return wrapped

            setattr(cls, name, make_wrapper(name, original))

        def measured_sorted(*args, **kw):
            start = time.perf_counter()
            result = builtins.sorted(*args, **kw)
            if sys._getframe(1).f_code.co_name == "enumerate":
                phases["result_sort"] = phases.get("result_sort", 0.) + time.perf_counter()-start
            return result

        enumerator.sorted = measured_sorted
    instance = cls(perms, graph["n_nodes"], **kwargs)
    gc.collect()
    if mode == "memory":
        import tracemalloc
        tracemalloc.start()
    started = time.perf_counter()
    result = instance.enumerate(labels)
    elapsed = time.perf_counter()-started
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_mb = rss / (1024 ** 2 if sys.platform == "darwin" else 1024)
    record = dict(backend=backend, case=[topology, c3, c4], mode=mode,
                  elapsed_s=elapsed, process_peak_rss_mb=rss_mb,
                  results=len(result), task_stats=instance.task_dedup_stats)
    if mode == "memory":
        current, peak = tracemalloc.get_traced_memory()
        snapshot = tracemalloc.take_snapshot()
        record["tracemalloc"] = dict(retained_bytes=current, peak_bytes=peak,
            top_lines=[dict(location=str(s.traceback[0]), bytes=s.size, allocations=s.count)
                       for s in snapshot.statistics("lineno")[:15]])
        tracemalloc.stop()
    if mode == "phase":
        del enumerator.sorted
        scheduler = phases["enumerate_subgroups"]
        search = phases["_search_compiled_task"]
        compile_time = phases["_compile_exclusions"]
        store = phases["_store_canonical_leaf"]
        lookup, insert = phases["_trie_lookup"], phases["_trie_insert"]
        sort, verify = phases.get("result_sort", 0.), phases["_verify_burnside"]
        record["instrumented_phases_s"] = dict(
            task_preparation_and_schedule=scheduler-search-compile_time,
            overlap_constraint_compilation=compile_time,
            search_excluding_leaf_storage=search-store,
            index_lookup=lookup, index_insert=insert,
            leaf_record_and_call_residual=store-lookup-insert,
            result_sort=sort, coverage_verification=verify,
            encode_decode_id_assignment_and_other=elapsed-scheduler-sort-verify)
        record["inclusive_method_times_s"] = phases
    if backend == "trie":
        stack = [instance._trie_root]
        count, size = 0, 0
        while stack:
            node = stack.pop()
            count += 1
            size += sys.getsizeof(node) + sys.getsizeof(node.children)
            stack.extend(node.children.values())
        index_stats = dict(nodes=count, index_container_bytes=size)
    else:
        index = instance._canonical_index
        index_stats = dict(keys=len(index), index_container_bytes=sys.getsizeof(index)
                           + sum(sys.getsizeof(key) for key in index))
    record["storage_structure"] = dict(**index_stats,
        result_container_bytes=sys.getsizeof(instance.classes) + sys.getsizeof(result)
        + sum(sys.getsizeof(r) + sys.getsizeof(r["coloring"]) for r in instance.classes.values()))
    record["output_sha256"] = hashlib.sha256(json.dumps(result, ensure_ascii=False,
                                      separators=(",", ":")).encode()).hexdigest()
    naive = math.prod(len(pools[d]) ** count for d, count in kwargs["topology_requirements"].items())
    assert sum(r["orbit_size"] for r in result) == naive
    assert all(instance._trie_lookup(tuple(r["coloring"])) is not None for r in result)
    return record


def launch(backend, case, mode):
    process = subprocess.run([sys.executable, __file__, "--worker", backend, *map(str, case), mode],
                             capture_output=True, text=True, check=True)
    return json.loads(process.stdout)


if __name__ == "__main__":
    if "--worker" in sys.argv:
        _, _, backend, topology, c3, c4, mode = sys.argv
        with contextlib.redirect_stdout(sys.stderr):
            record = worker(backend, topology, int(c3), int(c4), mode)
        print(json.dumps(record))
    else:
        records = []
        expected_digests = {}
        for case in CASES:
            results = {backend: [] for backend in BACKENDS}
            for repetition in range(3):
                for backend in BACKENDS[repetition:] + BACKENDS[:repetition]:
                    result = launch(backend, case, "timing")
                    key = tuple(case)
                    expected_digests.setdefault(key, result["output_sha256"])
                    assert result["output_sha256"] == expected_digests[key]
                    results[backend].append(result)
            records.append(dict(case=case, runs=results,
                medians_s={backend: statistics.median(r["elapsed_s"] for r in runs)
                           for backend, runs in results.items()},
                median_peak_rss_mb={backend: statistics.median(r["process_peak_rss_mb"] for r in runs)
                                    for backend, runs in results.items()}))
            print("Finished", case, records[-1]["medians_s"], file=sys.stderr)
        diagnostics = {}
        for mode in ["phase", "memory"]:
            diagnostics[mode] = []
            for backend in BACKENDS:
                result = launch(backend, CASES[-1], mode)
                assert result["output_sha256"] == expected_digests[CASES[-1]]
                diagnostics[mode].append(result)
                print("Finished", mode, backend, file=sys.stderr)
        print(json.dumps(dict(python=sys.version, platform=platform.platform(),
                              timings=records, diagnostics=diagnostics), indent=2))
