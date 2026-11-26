"""Same-store on/off comparison of completed-prefix subtree reuse.

Run after building the optional native core. Timings are five rotating warm
runs in one process, not directly comparable to the CASM fresh-process study.
Synthetic tasks expose a known overlap; real cases retain the usual schedule.
"""
import contextlib
import gc
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import statistics
import sys
import time

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import sympy.combinatorics  # Import cost is outside measurements.
from enumerator_native import HybridMOFEnumerator
from experiment._ours_worker import build_ligand_dict
from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
from symmetry import SymmetryEngine


class InterruptedSwapEnumerator(HybridMOFEnumerator):
    def enumerate_subgroups(self, labels):
        # Keep the real context/reset lifecycle; replace only the scheduler.
        from enumerator import MOFEnumerator
        from unittest.mock import patch
        def schedule(instance, atom_types):
            source = [[0, 1]] + [[i] for i in range(2, instance.n)]
            for partition, budget in ((source, 2 ** (instance.n - 2)),
                                      ([[i] for i in range(instance.n)], None)):
                instance._dfs_lattice_pruning(
                    0, partition, atom_types, [None] * instance.n, [1],
                    identity_start=len(instance.classes) if budget is not None else None,
                    identity_budget=budget)
        with patch.object(MOFEnumerator, "enumerate_subgroups", schedule):
            return super().enumerate_subgroups(labels)


class FirstTaskBudgetEnumerator(HybridMOFEnumerator):
    def _dfs_lattice_pruning(self, orbit_idx, *args, **kwargs):
        if orbit_idx == 0 and self.task_dedup_stats["tasks"] == 0:
            kwargs.update(identity_start=len(self.classes), identity_budget=3)
        return super()._dfs_lattice_pruning(orbit_idx, *args, **kwargs)


def compare(name, cls, perms, labels, options, naive):
    variants = [(backend, enabled) for backend in ("python", "cpp") for enabled in (False, True)]
    timings = {f"{b}_{int(e)}": [] for b, e in variants}
    reports = {}
    expected = None
    for repetition in range(6):  # First full rotation warms imports/allocators.
        order = variants[repetition % 4:] + variants[:repetition % 4]
        for backend, enabled in order:
            instance = cls(perms, len(perms[0]), backend=backend, prefix_trie=enabled,
                           verbose=False, max_total_results=None,
                           identity_orbit_limit=float("inf"), **options)
            gc.collect()
            with contextlib.redirect_stdout(io.StringIO()):
                started = time.perf_counter()
                records = instance.enumerate(labels)
                elapsed = time.perf_counter() - started
            assert instance.active_backend == backend
            if expected is None:
                expected = records
            assert records == expected, (name, backend, enabled)
            key = f"{backend}_{int(enabled)}"
            if repetition:
                timings[key].append(elapsed)
            reports[key] = dict(task_stats=instance.task_dedup_stats,
                                prefix_stats=instance.prefix_trie_stats,
                                native_kernel_s=sum(t["kernel_s"] for t in instance.native_task_reports))
    assert sum(r["orbit_size"] for r in expected) == naive
    return dict(case=name, sites=len(perms[0]), group_size=len(perms), naive=naive,
                records=len(expected), exact_records_order_and_metadata_match=True,
                sha256=hashlib.sha256(json.dumps(expected, separators=(",", ":")).encode()).hexdigest(),
                times_s=timings, median_s={key:statistics.median(t) for key,t in timings.items()},
                reports=reports)


def main():
    rows = []
    with contextlib.redirect_stdout(sys.stderr):
        for n in (12, 16, 18):
            perms = [tuple(range(n)), (1, 0) + tuple(range(2, n))]
            rows.append(compare(f"swap{n}_interrupted", InterruptedSwapEnumerator,
                                perms, ["A", "B"], {}, 2 ** n))
            print("Finished", rows[-1]["case"], rows[-1]["median_s"])
        for topology, c3, c4 in (("ffc", 3, 3), ("fjh", 1, 3)):
            graph = enrich_graph_data_with_topology(load_graph(str(ROOT/"data"/topology/f"{topology}_graph.json")))
            pools = build_ligand_dict(c3, c4)
            labels = sorted({label for pool in pools.values() for label in pool})
            degrees = [node["degree"] for node in graph["nodes"]]
            perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
            options = dict(topology_requirements=calculate_topology_requirements(graph),
                           user_selected_ligands=pools, node_degrees=degrees)
            import math
            naive = math.prod(len(pools[d]) for d in degrees)
            for kind, cls in (("exhaustive", HybridMOFEnumerator), ("first_budget3", FirstTaskBudgetEnumerator)):
                rows.append(compare(f"{topology}_{c3}_{c4}_{kind}", cls, perms, labels, options, naive))
                print("Finished", rows[-1]["case"], rows[-1]["median_s"])
    sources = ["search_prefix_trie.py", "enumerator_casm_inspired.py", "enumerator_native.py",
               "enumerator_storage.py", "enumerator_task_dedup.py", "enumerator.py",
               "native/enumeration_core.cpp", "native/libmof_enumeration.dylib",
               "validation/benchmark_prefix_trie_20261116.py"]
    # Resolve the Linux shared-library suffix when rerunning on that platform.
    sources = [name.replace(".dylib", ".so") if sys.platform != "darwin" else name for name in sources]
    output = dict(python=sys.version, platform=platform.platform(), native_abi=2,
                  repetitions=5, warmup_rotations=1, timing_scope="enumerate; same-process rotating warm runs; no RSS claim",
                  source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources},
                  rows=rows)
    path = ROOT/"validation/PREFIX_TRIE_BENCHMARK_2026-11-16.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(path)


if __name__ == "__main__":
    main()
