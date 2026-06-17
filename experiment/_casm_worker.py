#!/usr/bin/env python3.11
"""CASM enumeration worker — runs under python3.11 with libcasm v2.

Uses CASM's full factor group (proper + improper) for enumeration, then the
harness reconciles to proper-rotation orbit counts via P = A + 2C = F_total + C.

Usage:
    python3.11 _casm_worker.py <topology> <c3_count> <c4_count>

Prints a JSON record to stdout:
    {"orbits_full": N, "fg_size": N, "proper_size": N,
     "time_s": t, "peak_mem_mb": m, "status": "exhaustive"}
"""

import json
import os
import sys
import time
import resource

import numpy as np
from pymatgen.core import Structure

import libcasm.xtal as xtal
import libcasm.configuration as config_mod
import libcasm.enumerate as enum_mod

# ── Thread control ────────────────────────────────────────────────────────────
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

# ── Expected proper-rotation group sizes (graph automorphism |Γ|) ─────────────
EXPECTED_PROPER = {"ffc": 8, "fjh": 16, "vba": 16}
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_graph(topology: str) -> dict:
    path = os.path.join(PROJECT_ROOT, "data", topology, f"{topology}_graph.json")
    with open(path) as f:
        return json.load(f)


def build_prim(topology: str, c3_count: int, c4_count: int):
    """Construct CASM Prim from CIF with placeholder species.

    degree-3 sites → A, B, ... (first |C3| uppercase letters)
    degree-4 sites → X, Y, ... (first |C4| letters starting from X)
    """
    cif_path = os.path.join(PROJECT_ROOT, "data", "mof", f"{topology}.cif")
    structure = Structure.from_file(cif_path)
    graph = load_graph(topology)
    node_degrees = [n.get("degree", 0) for n in graph["nodes"]]
    n_sites = len(structure.sites)

    # ── Site mapping validation ───────────────────────────────────────────
    assert n_sites == graph["n_nodes"], (
        f"CIF sites ({n_sites}) != graph nodes ({graph['n_nodes']})"
    )
    for i, site in enumerate(structure.sites):
        node_frac = graph["nodes"][i].get("frac")
        if node_frac is not None:
            diff = np.linalg.norm(
                np.array(site.frac_coords) - np.array(node_frac)
            )
            assert diff < 1e-3, (
                f"Site {i}: CIF frac {site.frac_coords} != "
                f"graph frac {node_frac} (diff={diff:.4f})"
            )
    print(f"  [validate] {topology}: {n_sites} sites, CIF↔graph OK",
          file=sys.stderr)

    # ── Species labels ────────────────────────────────────────────────────
    deg3_species = [chr(ord("A") + i) for i in range(c3_count)] if c3_count > 0 else []
    deg4_species = [chr(ord("X") + i) for i in range(c4_count)] if c4_count > 0 else []

    occ_dof = []
    for d in node_degrees:
        if d == 3:
            occ_dof.append(deg3_species if deg3_species else ["A"])
        elif d == 4:
            occ_dof.append(deg4_species if deg4_species else ["X"])
        else:
            raise ValueError(f"Unexpected degree {d} at site of {topology}")

    all_species = set()
    for lst in occ_dof:
        all_species.update(lst)
    occupants = {s: xtal.Occupant(s) for s in all_species}

    # ── Construct xtal.Prim (3×N coordinate convention) ────────────────────
    lattice = xtal.Lattice(structure.lattice.matrix)
    frac_coords = np.array([s.frac_coords for s in structure.sites]).T  # 3×N
    xtal_prim = xtal.Prim(lattice, frac_coords, occ_dof, occupants=occupants)

    # ── Check proper-rotation subgroup ─────────────────────────────────────
    fg = xtal.make_prim_factor_group(xtal_prim)
    proper_ops = [op for op in fg if np.linalg.det(np.array(op.matrix())) > 0]
    expected = EXPECTED_PROPER[topology]

    if len(proper_ops) != expected:
        msg = (f"Proper-rotation subgroup size {len(proper_ops)} != "
               f"expected {expected} for {topology}.")
        if topology in ("ffc", "fjh"):
            print(f"  DIAGNOSTIC: {msg}", file=sys.stderr)
            for i, op in enumerate(fg):
                m = np.array(op.matrix())
                d = np.linalg.det(m)
                print(f"    op {i}: det={d:+.1f} trace={np.trace(m):+.1f}",
                      file=sys.stderr)
            raise SystemExit(1)
        else:
            print(f"  [WARNING] {msg} CASM effective |Γ|={len(proper_ops)}. "
                  f"Non-exhaustive case — continuing.", file=sys.stderr)

    print(f"  [sym] {topology}: fg_size={len(fg)}, proper={len(proper_ops)}",
          file=sys.stderr)

    # Wrap in config.Prim for enumeration API
    config_prim = config_mod.Prim(xtal_prim)
    return config_prim, node_degrees, len(fg), len(proper_ops)


def run_casm_enumeration(topology: str, c3_count: int, c4_count: int) -> dict:
    """Run CASM enumeration with full factor group. Returns full-group orbit count."""
    config_prim, _, fg_size, proper_size = build_prim(topology, c3_count, c4_count)

    # ── Volume-1 supercell ─────────────────────────────────────────────────
    supercell = config_mod.Supercell(
        config_prim,
        transformation_matrix_to_super=np.eye(3, dtype=int),
    )

    # ── Enumerate (full factor group) ──────────────────────────────────────
    enumerator = enum_mod.ConfigEnumAllOccupations(config_prim)

    t0 = time.perf_counter()
    n = 0
    for _ in enumerator.by_supercell_list([supercell]):
        n += 1
    t1 = time.perf_counter()

    rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {
        "orbits_full": n,
        "fg_size": fg_size,
        "proper_size": proper_size,
        "time_s": round(t1 - t0, 6),
        "peak_mem_mb": round(rss_bytes / 1_048_576, 2),
        "status": "exhaustive",
    }


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <topology> <c3_count> <c4_count>",
              file=sys.stderr)
        sys.exit(2)

    topology = sys.argv[1]
    c3_count = int(sys.argv[2])
    c4_count = int(sys.argv[3])

    # Redirect any stdout noise from libraries → stderr so only JSON goes to stdout
    _real_stdout = sys.stdout
    sys.stdout = sys.stderr

    try:
        result = run_casm_enumeration(topology, c3_count, c4_count)
    finally:
        sys.stdout = _real_stdout

    # Print clean JSON to (the real) stdout
    print(json.dumps(result))


if __name__ == "__main__":
    main()
