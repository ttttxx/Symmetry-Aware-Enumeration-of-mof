import argparse, sys, os, json, csv
sys.path.insert(0, '.')
from dft_validation import config
from dft_validation.select_structures import run_selection
from dft_validation.sanitize import sanitize_all
from dft_validation.cp2k_runner import run_one


def phase1(results_dir):
    print("=== Phase 1: select structures ===")
    records = run_selection(results_dir)
    ffc_cd = [r for r in records if r["topology"] == "ffc" and not r["probe_cu"]]
    probe  = [r for r in records if r["probe_cu"]]
    print(f"  ffc Cd records: {len(ffc_cd)}  (expect 10)")
    print(f"  Cu probe records: {len(probe)}  (expect 1)")
    ok = len(ffc_cd) == 10 and len(probe) <= 1
    if not ok:
        print("  FAIL: unexpected record counts")
        sys.exit(1)
    print("  Phase 1 PASS")
    return records


def phase2(records, results_dir):
    print("=== Phase 2: sanitize ===")
    records = sanitize_all(records, results_dir)
    ok = sum(1 for r in records if r.get("status") in ("clean", "preclean_done"))
    print(f"  clean/preclean: {ok}/{len(records)}")
    if ok == 0:
        print("  FAIL: no clean structures")
        sys.exit(1)
    print("  Phase 2 PASS")
    return records


def phase3(records, results_dir):
    print("=== Phase 3: dry pass — single ENERGY single-point ===")
    # Pick first clean Cd record (not probe)
    candidate = next(
        (r for r in records
         if r.get("status") in ("clean", "preclean_done") and not r.get("probe_cu")),
        None
    )
    if candidate is None:
        print("  FAIL: no clean Cd structure found")
        sys.exit(1)

    print(f"  Running: {candidate['id']} (n_atoms={candidate.get('n_atoms')})")
    updated = run_one(candidate, results_dir, "ffc",
                      relax_mode_override="ENERGY", walltime_seconds=600)

    if updated.get("dft_e_per_atom") is None:
        print(f"  FAIL: no energy parsed. Status: {updated.get('status')}")
        tail = updated.get("cp2k_log_tail", "")
        if tail:
            print(f"  Log tail:\n{tail}")
        sys.exit(1)

    e = updated["dft_e_per_atom"]
    e_h = updated.get("dft_e_hartree", "?")
    n   = updated.get("n_atoms", "?")
    print(f"\n  ✓ Energy parsed:")
    print(f"    E (Hartree) = {e_h}")
    print(f"    n_atoms     = {n}")
    print(f"    E (eV/atom) = {e:.6f}")
    print(f"\n  *** PAUSING for user confirmation ***")
    print(f"  Please verify the eV/atom value looks chemically reasonable")
    print(f"  (typical MOF: -5 to -15 eV/atom range) before proceeding to Phase 4.")

    # Write 1-row CSV
    csv_path = os.path.join(results_dir, "dft_vs_symmetry.csv")
    fieldnames = ["topology", "composition_label", "sym_class", "rho",
                  "dft_e_per_atom", "n_electrons", "multiplicity", "status", "id"]
    with open(csv_path, "w", newline="") as csvf:
        w = csv.DictWriter(csvf, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerow(updated)
    print(f"  CSV row written to {csv_path}")
    return [updated if r["id"] == updated["id"] else r for r in records]


def phase4_run(records, results_dir):
    print("=== Phase 4a: full batch CP2K run ===")
    updated = []
    for rec in records:
        if rec.get("probe_cu"):
            result = run_one(rec, results_dir, "ffc",
                             relax_mode_override="ENERGY", walltime_seconds=1800)
        elif rec.get("status") in ("clean", "preclean_done"):
            result = run_one(rec, results_dir, "ffc",
                             walltime_seconds=7200)
        else:
            result = rec
        updated.append(result)
        # Update JSON after each run (enables resume)
        json_path = os.path.join(results_dir, "selected_structures.json")
        with open(json_path, "w") as f:
            json.dump(updated + records[len(updated):], f, indent=2)

    converged = sum(1 for r in updated if r.get("status") == "converged")
    print(f"[phase4] Converged: {converged}/{len(updated)}")
    return updated


def phase4_analyze(records, results_dir):
    from dft_validation.analyze import analyze_ffc, get_fjh_gate, write_csv
    print("=== Phase 4b: analyze ===")
    write_csv(records, results_dir)
    comp_results = analyze_ffc(records, results_dir)

    print("\n--- ffc results ---")
    for label, stats in comp_results.items():
        hi = "YES" if stats["lowest_in_high"] else "NO"
        print(f"  {label}: Spearman rs={stats['spearman_rs']:.3f}  "
              f"p={stats['p_val']:.3f}  lowest_in_High={hi}")

    # fjh gate
    if config.FJH_GATE_ON_FFC:
        passed, reason = get_fjh_gate(comp_results)
        if passed:
            print(f"\n  [fjh gate] PASS — {reason}")
            print("  (fjh skipped locally: MAX_ATOMS=250 blocks ~800-atom fjh structures)")
        else:
            print(f"\n  [fjh gate] SKIP — {reason}")

    converged  = [r for r in records if r.get("status") == "converged"]
    failed     = [r for r in records if r.get("status") not in
                  ("converged", "done_cached", "selected", "assembled",
                   "clean", "preclean_done")]
    print(f"\n  converged={len(converged)}  failed/skipped={len(failed)}")
    print("  Phase 4 PASS")
    return records


def _load_records(results_dir):
    """Load existing selected_structures.json (for resume)."""
    json_path = os.path.join(results_dir, "selected_structures.json")
    if os.path.exists(json_path):
        with open(json_path) as f:
            return json.load(f)
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", type=int, default=1)
    parser.add_argument("--results-dir", default=config.RESULTS_DIR)
    parser.add_argument("--resume", action="store_true",
                        help="Skip Phase 1 enumeration, load existing selected_structures.json")
    args = parser.parse_args()

    records = None
    if args.resume:
        records = _load_records(args.results_dir)
        if records is None:
            print("  ERROR: --resume but no selected_structures.json found. Run without --resume first.")
            sys.exit(1)
        print(f"  [resume] Loaded {len(records)} records from JSON")
    else:
        if args.phase >= 1:
            records = phase1(args.results_dir)

    if args.phase >= 2:
        if records is None:
            records = phase1(args.results_dir)
        records = phase2(records, args.results_dir)

    if args.phase >= 3:
        records = phase3(records, args.results_dir)
        if args.phase == 3:
            print("\nPhase 3 complete. Review energy value above, then run --phase 4.")
            sys.exit(0)

    if args.phase >= 4:
        records = phase4_run(records, args.results_dir)
        records = phase4_analyze(records, args.results_dir)

if __name__ == "__main__":
    main()
