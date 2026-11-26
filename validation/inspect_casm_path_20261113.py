"""CASM 2.3.1 diagnostic: candidate count, proper subgroup, retained copies.

Run with python3.11; emits JSON and leaves original benchmark files unchanged.
Proper subgroup is only matched to ours after verifying exact site permutations.
"""
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import time

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiment._casm_worker import build_prim
import numpy as np
import libcasm.configuration as config
import libcasm.enumerate as enum
from libcasm.enumerate._enumerate import ConfigEnumAllOccupationsBase
from symmetry import SymmetryEngine
from main import load_graph, enrich_graph_data_with_topology


def check_case(topology, c3, c4):
    prim, degrees, fg_size, proper_size = build_prim(topology, c3, c4)
    supercell = config.Supercell(prim, transformation_matrix_to_super=np.eye(3, dtype=int))
    proper = []
    op, end = config.SupercellSymOp.begin(supercell), config.SupercellSymOp.end(supercell)
    while op != end:
        if np.linalg.det(op.to_symop().matrix()) > 0:
            proper.append(op.copy())
        op.next()
    graph = enrich_graph_data_with_topology(load_graph(str(ROOT / "data" / topology / f"{topology}_graph.json")))
    ours = set(SymmetryEngine(graph, verbose=False).compute_automorphism_group())
    # CASM pullback permutations are compared with inverses of ours' pushforward.
    ours_inverse = {tuple(p.index(i) for i in range(len(p))) for p in ours}
    casm_perms = {tuple(op.permute_index(i) for i in range(supercell.n_sites)) for op in proper}
    assert casm_perms == ours_inverse, "Proper group differs in its actual site action"
    n_sites = supercell.n_sites
    background = config.Configuration(supercell)
    naive = __import__("math").prod((c3 if d == 3 else c4) for d in degrees)

    # Exact existing public API; count-only and independent retained copies.
    public = []
    for store in (False, True):
        for _ in range(3):
            enumerator = enum.ConfigEnumAllOccupations(prim)
            saved = []
            count = 0
            start = time.perf_counter()
            for value in enumerator.by_supercell_list([supercell]):
                count += 1
                if store:
                    saved.append(value.copy())
            public.append(dict(store_copies=store, count=count, seconds=time.perf_counter()-start))

    # Diagnostic lower-level loop: same candidate generator, explicit subgroup,
    # no primitive check (volume 1); enumerate every full occupation once.
    restricted = []
    for _ in range(3):
        count, visited = 0, 0
        start = time.perf_counter()
        base = ConfigEnumAllOccupationsBase(background, set(range(n_sites)))
        while base.is_valid():
            visited += 1
            if config.is_canonical_configuration(base.value(), subgroup=proper):
                count += 1
            base.advance()
        restricted.append(dict(count=count, visited=visited, seconds=time.perf_counter()-start))
        assert visited == naive
    return dict(case=[topology, c3, c4], fg_size=fg_size, proper_size=proper_size,
                exact_permutation_group_match=True, naive=naive, public=public,
                proper_subgroup_lower_level=restricted)


if __name__ == "__main__":
    with contextlib.redirect_stdout(sys.stderr):
        records = [check_case(*case) for case in
                   [("ffc", 2, 2), ("ffc", 3, 3), ("fjh", 1, 2), ("fjh", 2, 1), ("fjh", 1, 3)]]
    print(json.dumps(records, indent=2))
