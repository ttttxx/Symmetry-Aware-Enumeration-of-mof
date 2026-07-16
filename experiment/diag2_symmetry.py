# !/usr/bin/env python3
"""Diagnostic 2: vba 8 vs 16 symmetry anomaly investigation.

Dumps vba's graph automorphism group Γ, runs SpacegroupAnalyzer at multiple
symprec on both the graph-baked structure and the raw CIF, and classifies
the discrepancy as Case A (tolerance artifact) or Case B (methodological).

Usage: python3 experiment/diag2_symmetry.py
"""

import json, os, sys
import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from symmetry import SymmetryEngine
from main import load_graph, enrich_graph_data_with_topology

def analyze_sga(structure, label, symprec_values=(0.01, 0.1, 0.3)):
    """Run SpacegroupAnalyzer at multiple symprec, report ops counts."""
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    print(f"  Formula: {structure.formula}, sites: {len(structure.sites)}")
    print(f"  Lattice: a={structure.lattice.a:.2f} b={structure.lattice.b:.2f} "
          f"c={structure.lattice.c:.2f}")

    for prec in symprec_values:
        try:
            sga = SpacegroupAnalyzer(structure, symprec=prec)
            ops = sga.get_symmetry_operations()
            total = len(ops)
            proper = sum(1 for op in ops if np.linalg.det(op.rotation_matrix) > 0)
            improper = total - proper
            sg = sga.get_space_group_symbol()
            sn = sga.get_space_group_number()
            print(f"  symprec={prec:4.2f}: total={total:3d}  proper={proper:3d}  "
                  f"improper={improper:3d}  spg={sg} (#{sn})")
        except Exception as e:
            print(f"  symprec={prec:4.2f}: ERROR {e}")
    print()

def dump_automorphism_group(graph_data):
    """Extract and dump the automorphism group Γ from SymmetryEngine."""
    print(f"\n{'='*60}")
    print(f"  Automorphism Group Γ — via SymmetryEngine")
    print(f"{'='*60}")

    sym = SymmetryEngine(graph_data, verbose=True)
    perms = sym.compute_automorphism_group()

    print(f"\n  |Γ| = {len(perms)}")
    print(f"  Node count: {sym.num_nodes}")

    # Build Structure for crystallographic reference
    lattice = sym.structure.lattice
    frac_coords = sym.structure.frac_coords

    # Get the crystallographic operations from SpacegroupAnalyzer
    sga_ops = []
    for prec in [0.01, 0.1, 0.3]:
        try:
            sga = SpacegroupAnalyzer(sym.structure, symprec=prec)
            sga_ops = sga.get_symmetry_operations()
            print(f"  [SGA symprec={prec}] structure ops: {len(sga_ops)} total, "
                  f"{sum(1 for o in sga_ops if np.linalg.det(o.rotation_matrix)>0)} proper")
            break  # use best
        except Exception:
            continue

    print(f"\n  Γ element analysis:")
    print(f"  {'idx':>3s}  {'det(R)':>7s}  {'type':>12s}  {'has (R,t)':>9s}  "
          f"  {'perm (first 8 of 32)':>30s}")
    print(f"  {'-'*3}  {'-'*7}  {'-'*12}  {'-'*9}  {'-'*30}")

    proper_count = 0
    improper_count = 0
    has_crystallographic = 0
    no_crystallographic = 0

    for idx, perm in enumerate(perms):
        # Try to find the crystallographic operation that maps to this permutation
        # The SymmetryEngine already mapped this — each perm came from an op with det>0.
        # But can we reconstruct the op?
        # SymmetryEngine code: for each op in sga_ops (already det>0 filtered):
        # → builds perm
        # So every perm in Γ was derived FROM a crystallographic op with det>0.
        # We need to find which one.
        # Strategy: try each op, apply to frac_coords, and see which gives this perm

        matched_op = None
        for op in sga_ops:
            if np.linalg.det(op.rotation_matrix) <= 0:
                continue
            candidate = []
            valid = True
            for i in range(sym.num_nodes):
                new_coord = op.operate(frac_coords[i])
                dists = lattice.get_all_distances(new_coord, frac_coords).flatten()
                match = np.where(dists < 0.8)[0]
                if len(match) == 0:
                    valid = False
                    break
                candidate.append(match[0])
            if valid and tuple(candidate) == tuple(perm):
                matched_op = op
                break

        # Determine properties
        rot_matrix = matched_op.rotation_matrix if matched_op else None
        det = np.linalg.det(rot_matrix) if rot_matrix is not None else None
        translation = matched_op.translation_vector if matched_op else None

        if det is not None:
            if det > 0:
                proper_count += 1
                op_type = "proper_rot"
            else:
                improper_count += 1
                op_type = "improper"
            has_crystallographic += 1
        else:
            op_type = "pure_graph"
            no_crystallographic += 1

        perm_str = str(list(perm)[:8]) if len(perm) > 8 else str(list(perm))
        det_str = f"{det:+.1f}" if det is not None else "N/A"
        has_op = "YES" if matched_op else "NO"
        rot_str = ""
        if rot_matrix is not None:
            rot_str = f"R={np.array2string(np.round(rot_matrix, 1), prefix=' ' * 30).replace(chr(10), '')}"
        trans_str = ""
        if translation is not None:
            trans_str = f" t=[{translation[0]:.3f},{translation[1]:.3f},{translation[2]:.3f}]"

        print(f"  {idx:3d}  {det_str:>7s}  {op_type:>12s}  {has_op:>9s}  "
              f"{perm_str:>30s}  {rot_str}{trans_str}")

    print(f"\n  Summary: |Γ|={len(perms)}  proper={proper_count}  "
          f"improper={improper_count}  has_crystallographic={has_crystallographic}  "
          f"pure_graph={no_crystallographic}")

    return perms, proper_count, improper_count, has_crystallographic, no_crystallographic

def analyze_cif_directly(cif_path, label):
    """Run SpacegroupAnalyzer directly on the CIF (what CASM sees)."""
    structure = Structure.from_file(cif_path)
    analyze_sga(structure, f"CIF: {label}")

def analyze_graph_structure(graph_data, label):
    """Run SpacegroupAnalyzer on the graph-baked structure (what SymmetryEngine sees)."""
    sym = SymmetryEngine(graph_data, verbose=False)
    # The structure is sym.structure
    analyze_sga(sym.structure, f"Graph-baked: {label}")

def main():
    vba_graph = load_graph(os.path.join(PROJECT_ROOT, "data", "vba", "vba_graph.json"))
    vba_graph = enrich_graph_data_with_topology(vba_graph)
    perms, proper, improper, has_crys, no_crys = dump_automorphism_group(vba_graph)

    analyze_cif_directly(
        os.path.join(PROJECT_ROOT, "data", "mof", "vba.cif"),
        "vba CIF (data/mof/vba.cif)"
    )

    vba_graph_raw = load_graph(os.path.join(PROJECT_ROOT, "data", "vba", "vba_graph.json"))
    analyze_graph_structure(vba_graph_raw, "vba graph-baked")

    for topo in ["ffc", "fjh"]:
        analyze_cif_directly(
            os.path.join(PROJECT_ROOT, "data", "mof", f"{topo}.cif"),
            f"{topo} CIF (data/mof/{topo}.cif)"
        )

    print(f"\n{'='*60}")
    print(f"  CLASSIFICATION")
    print(f"{'='*60}")

    if improper > 0 or no_crys > 0:
        print(f"  → Case B (methodological issue):")
        print(f"    Γ contains {improper} improper rotation(s) or {no_crys} pure-graph"
              f" automorphism(s) without crystallographic (R,t) correspondence.")
        print(f"    The graph automorphism group differs from the crystallographic symmetry.")
    elif has_crys == len(perms):
        print(f"  → Case A (tolerance artifact):")
        print(f"    All {len(perms)} Γ elements have proper-rotation (R,t) correspondences.")
        print(f"    Discrepancy in CASM is due to tolerance/symprec differences or")
        print(f"    coordinate precision between CIF and graph-baked structure.")
    else:
        print(f"  → Mixed: {has_crys} with (R,t), {no_crys} without. Need further investigation.")

if __name__ == "__main__":
    main()
