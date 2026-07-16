# !/usr/bin/env python3.11
"""Diagnostic 3: CASM throughput measurement and vba extrapolation.

Measures CASM enumeration throughput on the largest feasible case (fjh C3=1,C4=3)
and extrapolates vba wall-clock lower bounds.
"""
import json, os, sys, time
import numpy as np
from pymatgen.core import Structure

import libcasm.xtal as xtal
import libcasm.configuration as config_mod
import libcasm.enumerate as enum_mod

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def measure_throughput(topology, c3, c4):
    """Measure CASM enumeration throughput (naive_configs/sec and orbits/sec)."""
    cif_path = os.path.join(PROJECT_ROOT, "data", "mof", f"{topology}.cif")
    structure = Structure.from_file(cif_path)
    graph = json.load(open(
        os.path.join(PROJECT_ROOT, "data", topology, f"{topology}_graph.json")))
    nd = [n.get("degree", 0) for n in graph["nodes"]]

    deg3_species = [chr(ord("A") + i) for i in range(c3)] if c3 else ["A"]
    deg4_species = [chr(ord("X") + i) for i in range(c4)] if c4 else ["X"]
    occ_dof = [deg3_species if d == 3 else deg4_species for d in nd]
    all_sp = set()
    for lst in occ_dof:
        all_sp.update(lst)

    lattice = xtal.Lattice(structure.lattice.matrix)
    frac = np.array([s.frac_coords for s in structure.sites]).T
    xtal_prim = xtal.Prim(lattice, frac, occ_dof,
                          occupants={s: xtal.Occupant(s) for s in all_sp})
    config_prim = config_mod.Prim(xtal_prim)
    supercell = config_mod.Supercell(config_prim,
                                      transformation_matrix_to_super=np.eye(3, dtype=int))

    # Compute naive space
    naive = 1
    for d in nd:
        n_choices = len(deg3_species) if d == 3 else len(deg4_species)
        naive *= n_choices

    fg = xtal.make_prim_factor_group(xtal_prim)
    fg_size = len(fg)

    # Time the enumeration
    e = enum_mod.ConfigEnumAllOccupations(config_prim)
    t0 = time.perf_counter()
    n = 0
    for _ in e.by_supercell_list([supercell]):
        n += 1
    t1 = time.perf_counter()
    elapsed = t1 - t0

    throughput_naive = naive / elapsed if elapsed > 0 else float("inf")
    throughput_orbits = n / elapsed if elapsed > 0 else float("inf")

    return {
        "topology": topology,
        "c3": c3, "c4": c4,
        "n_sites": len(nd),
        "naive_space": naive,
        "fg_size": fg_size,
        "orbits_full": n,
        "time_s": elapsed,
        "naive_per_sec": throughput_naive,
        "orbits_per_sec": throughput_orbits,
    }

def extrapolate_vba(throughput_naive_per_sec):
    """Extrapolate vba wall-clock from naive-space throughput."""
    results = []
    for c4, label in [(2, "vba C4=2"), (3, "vba C4=3")]:
        naive = c4 ** 32  # 32 nodes all deg-4
        # Use CASM's |Γ|=8 for estimate
        orbits_est = naive / 8
        time_est_sec = naive / throughput_naive_per_sec if throughput_naive_per_sec > 0 else float("inf")
        time_est_hours = time_est_sec / 3600
        time_est_years = time_est_hours / 8760

        results.append({
            "case": label,
            "naive_space": naive,
            "naive_log10": np.log10(naive),
            "orbits_est": orbits_est,
            "time_est_sec": time_est_sec,
            "time_est_hours": time_est_hours,
            "time_est_years": time_est_years,
            "intractable_by_30min": time_est_sec > 1800,
            "intractable_by_8GB": True,  # Memory also a concern
        })
    return results

def main():
    print("=" * 60)
    print("  CASM Throughput Measurement")
    print("=" * 60)

    # Measure on the largest feasible case
    print("\n--- fjh C3=1,C4=3 (most comparable to vba) ---")
    r = measure_throughput("fjh", 1, 3)
    for k, v in r.items():
        print(f"  {k}: {v}")

    print(f"\n  Throughput: {r['naive_per_sec']:.1f} naive configs/sec")
    print(f"  Throughput: {r['orbits_per_sec']:.1f} orbit reps/sec")

    # Also measure on a smaller case for comparison
    print("\n--- fjh C3=1,C4=2 (smaller, for comparison) ---")
    r2 = measure_throughput("fjh", 1, 2)
    print(f"  Throughput: {r2['naive_per_sec']:.1f} naive configs/sec")
    print(f"  Throughput: {r2['orbits_per_sec']:.1f} orbit reps/sec")

    # Use the more conservative (slower) throughput for extrapolation
    tp = min(r["naive_per_sec"], r2["naive_per_sec"])
    print(f"\n  Conservative throughput: {tp:.1f} naive configs/sec")

    # Extrapolate
    print(f"\n{'='*60}")
    print(f"  vba Extrapolation (throughput={tp:.1f} naive/sec, |Γ_CASM|=8)")
    print(f"{'='*60}")

    for c4 in [2, 3]:
        naive = c4 ** 32
        orbits_est = naive / 8
        time_est = naive / tp
        print(f"\n  vba C4={c4}:")
        print(f"    naive space = {c4}^32 = {naive:.2e}")
        print(f"    est. orbits (|Γ|=8) = {orbits_est:.2e}")
        print(f"    est. wall time = {time_est:.1f}s = {time_est/3600:.1f}h "
              f"= {time_est/8760:.1f} years")
        print(f"    intractable? {'YES' if time_est > 1800 else 'MAYBE'} "
              f"(threshold: 30min = 1800s)")

    # Check: is the 30min watchdog appropriate or too short?
    print(f"\n{'='*60}")
    print(f"  Watchdog Assessment")
    print(f"{'='*60}")
    print(f"  Current harness watchdog: 30 min (1800s) or 8 GB RSS")
    c4_2_time = (2**32) / tp
    print(f"  vba C4=2 would need ~{c4_2_time/3600:.1f} hours")
    print(f"  → Would trigger 30min watchdog: {'YES' if c4_2_time > 1800 else 'NO'}")
    if c4_2_time > 1800:
        print(f"  → Conclusion: vba C4=2 is genuinely intractable, not post-hoc cutoff")
    c4_3_time = (3**32) / tp
    print(f"  vba C4=3 would need ~{c4_3_time/8760:.1f} years")
    print(f"  → Clearly intractable by many orders of magnitude")

if __name__ == "__main__":
    main()
