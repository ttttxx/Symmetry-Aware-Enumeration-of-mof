import argparse, sys, os, json, re
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
sys.path.insert(0, '.')
from dft_validation import config

def _parse_cp2k_output(out_path: str, n_atoms: int) -> tuple[float | None, float | None]:
    if not os.path.exists(out_path):
        return None, None
    with open(out_path, errors="replace") as f:
        text = f.read()
    matches = re.findall(
        r"ENERGY\| Total FORCE_EVAL[^\n]*\s+(-?\d+\.\d+)", text
    )
    if not matches:
        return None, None
    e_hartree = float(matches[-1])
    e_per_atom = e_hartree * config.HARTREE_TO_EV / max(1, n_atoms)
    return e_per_atom, e_hartree

def _walltime_for(topology: str, default: int = 7200) -> int:
    wt = getattr(config, "WALLTIME_SECONDS", None)
    if isinstance(wt, dict):
        return wt.get(topology, default)
    if isinstance(wt, int):
        return wt
    return default

def _workers_for(topology: str) -> int:
    if not getattr(config, "DFT_PARALLEL", False):
        return 1
    workers = getattr(config, "DFT_MAX_WORKERS", 1)
    if isinstance(workers, dict):
        workers = workers.get(topology, 1)
    try:
        workers = int(workers)
    except (TypeError, ValueError):
        workers = 1
    if workers < 1:
        return 1
    omp = max(1, int(getattr(config, "CP2K_OMP", 1)))
    cpu_cap = max(1, (os.cpu_count() or omp) // omp)
    return max(1, min(workers, cpu_cap))

def _mem_for(topology: str, workers: int) -> float:
    per_job = getattr(config, "DFT_JOB_MEM_GB", None)
    if isinstance(per_job, dict) and topology in per_job:
        return float(per_job[topology])
    if isinstance(per_job, (int, float)):
        return float(per_job)
    return float(getattr(config, "MAX_MEM_GB", 100)) / max(1, workers)

def _eligible_for_dft(rec: dict) -> bool:
    return rec.get("probe_cu") or rec.get("status") in ("clean", "preclean_done")

def _run_record(rec: dict, results_dir: str, workers: int) -> dict:
    from dft_validation.cp2k_runner import run_one

    topo = rec.get("topology", "ffc")
    kwargs = dict(
        walltime_seconds=_walltime_for(topo),
        max_mem_gb=_mem_for(topo, workers),
    )
    if rec.get("probe_cu"):
        kwargs["relax_mode_override"] = "ENERGY"
    return run_one(rec, results_dir, topo, **kwargs)

def harvest_existing_outputs(records: list, results_dir: str) -> list:
    """Parse existing cp2k_run.out files into records without launching CP2K."""
    updated = []
    harvested = 0
    missing = 0
    for rec in records:
        topo = rec.get("topology", "ffc")
        rec_id = rec["id"]
        out_path = os.path.join(results_dir, topo, "runs", rec_id, "cp2k_run.out")
        e_per_atom, e_hartree = _parse_cp2k_output(out_path, rec.get("n_atoms") or 1)
        if e_per_atom is None:
            missing += 1
            updated.append(rec)
            continue
        harvested += 1
        relaxed_cif = os.path.join(results_dir, topo, "relaxed", f"{rec_id}.cif")
        new_rec = {
            **rec,
            "dft_e_per_atom": e_per_atom,
            "dft_e_hartree": e_hartree,
            "status": "converged",
        }
        if os.path.exists(relaxed_cif):
            new_rec["relaxed_cif"] = relaxed_cif
        updated.append(new_rec)

    json_path = os.path.join(results_dir, "selected_structures.json")
    with open(json_path, "w") as f:
        json.dump(updated, f, indent=2)
    print(f"  [harvest] parsed energies from existing CP2K outputs: "
          f"{harvested}/{len(records)}; missing/no-energy={missing}")
    print(f"  [harvest] updated {json_path}")
    return updated

def _merge_records_preserve_existing(existing: list, proposed: list) -> list:
    existing_by_id = {r["id"]: r for r in existing}
    proposed_ids = {r["id"] for r in proposed}

    merged = []
    kept = 0
    added = 0
    for rec in proposed:
        old = existing_by_id.get(rec["id"])
        if old is not None:
            merged.append(old)
            kept += 1
        else:
            merged.append(rec)
            added += 1

    orphaned = [r for r in existing if r["id"] not in proposed_ids]
    merged.extend(orphaned)
    if orphaned:
        print(f"  [expand] kept {len(orphaned)} existing records not in new selection")
    print(f"  [expand] preserved existing={kept}, added new={added}, total={len(merged)}")
    return merged

def expand_selection_preserving_results(existing: list, results_dir: str) -> list:
    """Re-run selection with current config, merge by id, then sanitize new records."""
    from dft_validation.select_structures import run_selection
    from dft_validation.sanitize import sanitize_all

    print("=== Expand selection: preserve existing DFT results, add new candidates ===")
    proposed = run_selection(results_dir)
    merged = _merge_records_preserve_existing(existing, proposed)
    merged = sanitize_all(merged, results_dir)
    json_path = os.path.join(results_dir, "selected_structures.json")
    with open(json_path, "w") as f:
        json.dump(merged, f, indent=2)
    print(f"  [expand] wrote merged selection to {json_path}")
    return merged

def phase1(results_dir):
    from dft_validation.select_structures import run_selection

    print("=== Phase 1: select structures ===")
    records = run_selection(results_dir)
    assembled = [r for r in records if r.get("status") == "assembled"]
    by_topo = Counter(r["topology"] for r in assembled)
    print(f"  assembled by topology: {dict(by_topo)}")
    if not assembled:
        print("  FAIL: nothing assembled")
        sys.exit(1)
    # Each configured topology should contribute at least one assembled structure.
    for t in config.VALIDATION_TARGETS:
        if by_topo.get(t["topology"], 0) == 0:
            print(f"  WARNING: topology {t['topology']} produced 0 assembled structures "
                  f"(check distinct-rho / MIN_RHO settings)")
    print("  Phase 1 PASS")
    return records

def phase2(records, results_dir):
    from dft_validation.sanitize import sanitize_all

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
    """Local smoke test: 1 ENERGY single-point at SMOKE_CUTOFF_RY (non-production)."""
    from dft_validation.cp2k_runner import run_one

    print("=== Phase 3: smoke — 1 ENERGY single-point ===")
    candidate = next(
        (r for r in records
         if r.get("status") in ("clean", "preclean_done") and not r.get("probe_cu")),
        None
    )
    if candidate is None:
        print("  FAIL: no clean structure found")
        sys.exit(1)

    topo = candidate["topology"]
    print(f"  Target:  {candidate['id']} (topology={topo}, n_atoms={candidate.get('n_atoms')})")
    print(f"  Cutoff:  {config.SMOKE_CUTOFF_RY} Ry  [non-production]")
    updated = run_one(candidate, results_dir, topo,
                      relax_mode_override="ENERGY",
                      walltime_seconds=3600,
                      cutoff_ry=config.SMOKE_CUTOFF_RY)

    if updated.get("dft_e_per_atom") is None:
        print(f"  FAIL: no energy parsed. Status: {updated.get('status')}")
        tail = updated.get("cp2k_log_tail", "")
        if tail:
            print(f"  Log tail:\n{tail}")
        sys.exit(1)

    print(f"\n  ✓ E (Hartree)={updated.get('dft_e_hartree')}  "
          f"n_atoms={updated.get('n_atoms')}  E(eV/atom)={updated['dft_e_per_atom']:.6f}")
    print("  Smoke test PASS")
    return [updated if r["id"] == updated["id"] else r for r in records]

def phase4_run(records, results_dir, only_id=None):
    """Full batch CP2K run. Uses each record's OWN topology (no hard-coding).
    ENERGY single-points; per-topology walltime; resumable JSON after each run.
    """
    print("=== Phase 4a: full batch CP2K run ===")
    if not config.RUN_DFT_LOCALLY:
        print("  [SKIP] RUN_DFT_LOCALLY=False.")
        return records

    if only_id:
        if not any(r["id"] == only_id for r in records):
            print(f"  ERROR: --only {only_id!r} not found.")
            sys.exit(1)
        print(f"  [--only] Running single structure: {only_id}")

    all_records = list(records)
    json_path = os.path.join(results_dir, "selected_structures.json")

    if only_id:
        for i, rec in enumerate(all_records):
            if rec["id"] != only_id:
                continue
            all_records[i] = _run_record(rec, results_dir, workers=1)
            with open(json_path, "w") as f:
                json.dump(all_records, f, indent=2)
            break
    else:
        indices_by_topo = defaultdict(list)
        for i, rec in enumerate(all_records):
            if not _eligible_for_dft(rec):
                print(f"  [resume] {rec['id']}: skip status={rec.get('status')}")
                continue
            indices_by_topo[rec.get("topology", "ffc")].append(i)

        topo_order = [t["topology"] for t in config.VALIDATION_TARGETS]
        topo_order.extend(t for t in indices_by_topo if t not in topo_order)

        for topo in topo_order:
            indices = indices_by_topo.get(topo, [])
            if not indices:
                continue
            workers = _workers_for(topo)
            mem = _mem_for(topo, workers)
            active_threads = workers * int(getattr(config, "CP2K_OMP", 1))
            print(f"  [batch:{topo}] {len(indices)} structures; "
                  f"workers={workers}; OMP/job={config.CP2K_OMP}; "
                  f"active_threads≈{active_threads}; mem/job≤{mem:g} GB")

            if workers == 1:
                for i in indices:
                    all_records[i] = _run_record(all_records[i], results_dir, workers)
                    with open(json_path, "w") as f:
                        json.dump(all_records, f, indent=2)
                continue

            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(_run_record, all_records[i], results_dir, workers): i
                    for i in indices
                }
                for future in as_completed(futures):
                    i = futures[future]
                    try:
                        result = future.result()
                    except Exception as exc:
                        result = {**all_records[i], "status": f"runner_exception: {exc}"}
                    all_records[i] = result
                    print(f"  [batch:{topo}] finished {result['id']} "
                          f"status={result.get('status')}")
                    with open(json_path, "w") as f:
                        json.dump(all_records, f, indent=2)

    converged   = sum(1 for r in all_records if r.get("status") == "converged")
    done_cached = sum(1 for r in all_records if r.get("status") == "done_cached")
    print(f"  converged={converged}  done_cached={done_cached}  total={len(all_records)}")
    return all_records

def phase4_analyze(records, results_dir):
    """CSV + per-topology rho–E figures + per-composition Spearman for ALL topologies."""
    from dft_validation.analyze import analyze_all, get_fjh_gate, write_csv
    print("=== Phase 4b: analyze ===")
    write_csv(records, results_dir)
    by_topo = analyze_all(records, results_dir)

    for topo, comp_results in by_topo.items():
        print(f"\n--- {topo} results ---")
        for label, stats in comp_results.items():
            rs = stats["spearman_rs"]
            rs_str = f"{rs:.3f}" if rs == rs else "n/a"
            p_str  = f"{stats['p_val']:.3f}" if stats['p_val'] == stats['p_val'] else "n/a"
            hi = "YES" if stats["lowest_at_max_rho"] else "NO"
            print(f"  {label}: n={stats['n']}  Spearman rs={rs_str}  p={p_str}  "
                  f"lowest_at_max_rho={hi}")
        passed, reason = get_fjh_gate(comp_results)
        if passed is None:
            verdict = "PENDING"
        else:
            verdict = "SUPPORTS hypothesis" if passed else "does NOT support"
        print(f"  [{topo} qualitative] {verdict} — {reason}")

    failed = [r for r in records if r.get("status") not in
              ("converged", "done_cached", "selected", "assembled",
               "clean", "preclean_done")]
    if failed:
        print(f"\n  failed/anomalous statuses: "
              f"{Counter(r.get('status') for r in failed)}")
    print("  Phase 4 PASS")
    return records

def print_energy_summary(records: list) -> None:
    """Compact table for checking whether current energies look sane."""
    rows = [r for r in records if r.get("dft_e_per_atom") is not None]
    if not rows:
        print("\n=== Energy summary ===")
        print("  no parsed DFT energies yet")
        return

    print("\n=== Energy summary ===")
    for topo in sorted({r.get("topology", "ffc") for r in rows}):
        topo_rows = sorted(
            [r for r in rows if r.get("topology", "ffc") == topo],
            key=lambda r: (r.get("composition_label", ""), -r.get("rho", 0), r["id"]),
        )
        vals = [r["dft_e_per_atom"] for r in topo_rows]
        print(f"  {topo}: n={len(topo_rows)}  "
              f"E_range=[{min(vals):.6f}, {max(vals):.6f}] eV/atom")
        for r in topo_rows:
            print(f"    {r['id']:20s} {r.get('composition_label',''):8s} "
                  f"rho={r.get('rho', float('nan')):7.4f} "
                  f"{r.get('sym_class',''):4s} "
                  f"E={r['dft_e_per_atom']: .6f} eV/atom")

def _load_records(results_dir):
    json_path = os.path.join(results_dir, "selected_structures.json")
    if os.path.exists(json_path):
        with open(json_path) as f:
            return json.load(f)
    return None

def main():
    parser = argparse.ArgumentParser(
        description="DFT validation pipeline: select→sanitize→smoke→DFT→analyze"
    )
    parser.add_argument("--phase", type=int, default=1, help="Run up to this phase (1-4)")
    parser.add_argument("--results-dir", default=config.RESULTS_DIR)
    parser.add_argument("--resume", action="store_true",
                        help="Skip Phase 1 enumeration; load existing selected_structures.json")
    parser.add_argument("--only", metavar="ID", help="Phase 4: run only this structure id")
    parser.add_argument("--analyze-only", action="store_true",
                        help="Skip Phase 4a run; only produce CSV+figures from existing JSON")
    parser.add_argument("--harvest-only", action="store_true",
                        help="Parse existing cp2k_run.out files, update JSON, then analyze")
    parser.add_argument("--expand-selection", action="store_true",
                        help="Re-run selection with current config, preserve existing results, add new records")
    parser.add_argument("--summary", action="store_true",
                        help="Print compact parsed-energy table after analyze/harvest")
    args = parser.parse_args()

    skip_to_phase4 = (config.START_FROM == "clean" and args.resume)

    records = None
    if args.resume or skip_to_phase4:
        records = _load_records(args.results_dir)
        if records is None:
            print("  ERROR: --resume but no selected_structures.json found.")
            sys.exit(1)
        print(f"  [resume] Loaded {len(records)} records from JSON")

    if not skip_to_phase4:
        if args.phase >= 1 and records is None:
            records = phase1(args.results_dir)
        if args.phase >= 2:
            if records is None:
                records = phase1(args.results_dir)
            records = phase2(records, args.results_dir)
        if args.phase >= 3:
            records = phase3(records, args.results_dir)
            if args.phase == 3:
                print("\nPhase 3 complete. Review energy above, then run --phase 4.")
                sys.exit(0)

    if args.phase >= 4 or skip_to_phase4:
        if args.expand_selection:
            if records is None:
                records = _load_records(args.results_dir)
                if records is None:
                    print("  ERROR: --expand-selection requires existing selected_structures.json")
                    sys.exit(1)
            records = harvest_existing_outputs(records, args.results_dir)
            records = expand_selection_preserving_results(records, args.results_dir)
        elif args.harvest_only:
            if records is None:
                records = _load_records(args.results_dir)
                if records is None:
                    print("  ERROR: --harvest-only requires existing selected_structures.json")
                    sys.exit(1)
            records = harvest_existing_outputs(records, args.results_dir)
        elif args.analyze_only:
            if records is None:
                records = _load_records(args.results_dir)
                if records is None:
                    print("  ERROR: --analyze-only requires existing selected_structures.json")
                    sys.exit(1)
        else:
            records = phase4_run(records, args.results_dir, only_id=args.only)
        records = phase4_analyze(records, args.results_dir)
        if args.summary:
            print_energy_summary(records)

if __name__ == "__main__":
    main()
