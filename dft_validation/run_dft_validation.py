import argparse, sys, os
sys.path.insert(0, '.')
from dft_validation import config
from dft_validation.select_structures import run_selection


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", type=int, default=1)
    parser.add_argument("--results-dir", default=config.RESULTS_DIR)
    args = parser.parse_args()

    if args.phase >= 1:
        records = phase1(args.results_dir)
    # Phases 2–4 added in later tasks

if __name__ == "__main__":
    main()
