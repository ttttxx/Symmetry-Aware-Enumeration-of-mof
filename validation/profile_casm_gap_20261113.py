"""Read-only benchmark diagnostics; does not update paper or original CSV.

Run with the project Python. Each case is a fresh process; cold baseline is
followed by three warm baselines and three identity-only diagnostic runs.
Identity-only preserves the production DFS and leaf processing, but omits
subgroup scheduling. It is an exhaustive-mode ablation, not a replacement for
the symmetry-prioritized bounded-output algorithm.
"""
import contextlib
import cProfile
import io
import json
import os
from pathlib import Path
import pstats
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
CASES = [("ffc", 2, 2), ("ffc", 3, 3), ("fjh", 1, 2), ("fjh", 2, 1), ("fjh", 1, 3)]


def run_case(topology, c3, c4):
    from experiment._ours_worker import build_ligand_dict
    from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
    from symmetry import SymmetryEngine
    from enumerator import MOFEnumerator

    graph = enrich_graph_data_with_topology(load_graph(str(ROOT / "data" / topology / f"{topology}_graph.json")))
    degrees = [node["degree"] for node in graph["nodes"]]
    requirements = calculate_topology_requirements(graph)
    perms = SymmetryEngine(graph, verbose=False).compute_automorphism_group()
    pools = build_ligand_dict(c3, c4)
    labels = sorted({label for pool in pools.values() for label in pool})

    class IdentityOnly(MOFEnumerator):
        def enumerate_subgroups(self, atom_types):
            self._last_atom_types = atom_types
            self._dfs_lattice_pruning(0, [[i] for i in range(self.n)], atom_types,
                                      [None] * self.n, list(range(1, len(self.perms))))

    class Instrumented(MOFEnumerator):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.tasks = []
            self.leaves = 0
            self.method_seconds = {}

        def _dfs_lattice_pruning(self, orbit_idx, *args, **kwargs):
            if orbit_idx != 0:
                return super()._dfs_lattice_pruning(orbit_idx, *args, **kwargs)
            start = time.perf_counter()
            if not self.tasks:
                self.setup_s = start - self.start
            before, leaves = len(self.classes), self.leaves
            result = super()._dfs_lattice_pruning(orbit_idx, *args, **kwargs)
            self.tasks.append(dict(orbits=len(args[0]), seconds=time.perf_counter()-start,
                                   leaves=self.leaves-leaves, new=len(self.classes)-before))
            return result

        def _try_add_result(self, coloring):
            self.leaves += 1
            return super()._try_add_result(coloring)

    def timed_method(name):
        original = getattr(MOFEnumerator, name)
        def wrapped(self, *args, **kwargs):
            start = time.perf_counter()
            result = original(self, *args, **kwargs)
            self.method_seconds[name] = self.method_seconds.get(name, 0) + time.perf_counter()-start
            return result
        return wrapped

    for name in ("_get_canonical", "_compute_stabilizer_and_orbit", "_trie_lookup", "_trie_insert"):
        setattr(Instrumented, name, timed_method(name))

    def make(cls):
        return cls(perms, graph["n_nodes"], verbose=False,
                   topology_requirements=requirements, user_selected_ligands=pools,
                   node_degrees=degrees, max_total_results=None, identity_orbit_limit=float("inf"))

    expected = None
    def execute(cls):
        nonlocal expected
        obj = make(cls)
        start = time.perf_counter()
        result = obj.enumerate(labels)
        elapsed = time.perf_counter() - start
        signature = {tuple(r["coloring"]): (r["stabilizer_size"], r["orbit_size"]) for r in result}
        if expected is None:
            expected = signature
        assert signature == expected, "Ablation changed representatives or orbit metadata"
        assert sum(r["orbit_size"] for r in result) == __import__("math").prod(len(pools[d]) ** n for d, n in requirements.items())
        return elapsed

    cold = execute(MOFEnumerator)
    warm = [execute(MOFEnumerator) for _ in range(3)]
    identity = [execute(IdentityOnly) for _ in range(3)]
    obj = make(Instrumented)
    obj.start = time.perf_counter()
    obj.enumerate(labels)
    profile = cProfile.Profile()
    profile.enable()
    execute(MOFEnumerator)
    profile.disable()
    buffer = io.StringIO()
    pstats.Stats(profile, stream=buffer).strip_dirs().sort_stats("cumulative").print_stats(20)
    return dict(case=[topology, c3, c4], group_size=len(perms), orbits=len(expected),
                cold_s=cold, warm_s=warm, identity_only_s=identity,
                instrumented_setup_s=obj.setup_s, instrumented_tasks=obj.tasks,
                instrumented_method_s=obj.method_seconds, profile=buffer.getvalue())


if __name__ == "__main__":
    if len(sys.argv) == 4:
        with contextlib.redirect_stdout(sys.stderr):
            record = run_case(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
        print(json.dumps(record))
    else:
        results = []
        for topology, c3, c4 in CASES:
            process = subprocess.run([sys.executable, __file__, topology, str(c3), str(c4)],
                                     capture_output=True, text=True, check=True)
            results.append(json.loads(process.stdout))
            print(f"Finished {topology} {c3}/{c4}", file=sys.stderr)
        print(json.dumps(results, indent=2))
