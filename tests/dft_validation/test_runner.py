import sys, os, re
sys.path.insert(0, '.')
import pytest


def test_input_contains_pbe():
    import glob
    cifs = glob.glob("results/dft_validation/ffc/clean/*.cif")
    if not cifs:
        pytest.skip("Phase 2 outputs not present")
    import json
    recs = json.load(open("results/dft_validation/selected_structures.json"))
    rec = next(r for r in recs if r.get("status") in ("clean", "preclean_done", "converged", "done_cached")
               and not r.get("probe_cu"))
    from dft_validation.cp2k_runner import generate_cp2k_input
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        inp = generate_cp2k_input(rec, td, "ffc", relax_mode_override="ENERGY")
        text = open(inp).read()
    assert "PBE" in text
    assert "DFTD3" in text
    assert "ENERGY" in text
    assert "DZVP-MOLOPT-PBE-GTH" in text


def test_input_charge_zero():
    import glob, json, tempfile
    cifs = glob.glob("results/dft_validation/ffc/clean/*.cif")
    if not cifs:
        pytest.skip("Phase 2 outputs not present")
    recs = json.load(open("results/dft_validation/selected_structures.json"))
    rec = next(r for r in recs if r.get("status") in ("clean", "preclean_done", "converged", "done_cached"))
    from dft_validation.cp2k_runner import generate_cp2k_input
    with tempfile.TemporaryDirectory() as td:
        inp = generate_cp2k_input(rec, td, "ffc", relax_mode_override="ENERGY")
        text = open(inp).read()
    assert re.search(r"CHARGE\s+0", text)


def test_energy_parse():
    from dft_validation.cp2k_runner import parse_energy
    fake_out = """
 SCF run converged in    15 steps
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:            -5812.123456789012
"""
    e = parse_energy(fake_out, n_atoms=200)
    assert e is not None
    expected = -5812.123456789012 * 27.211386 / 200
    assert abs(e - expected) < 1e-6


def test_relaxation_completion_distinguishes_partial_restart():
    from dft_validation.cp2k_runner import parse_relaxation_completion
    partial_out = """
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:            -1032.6865048233376
 --------------------------
 OPTIMIZATION STEP:     41
 --------------------------
 OPT| Maximum step size is converged                                          NO
 OPT| RMS step size is converged                                              NO
 OPT| Maximum gradient is converged                                           NO
 OPT| RMS gradient is converged                                               NO
 OPT| Pressure is converged                                                   NO
"""
    partial = parse_relaxation_completion(partial_out, "CELL_OPT")
    assert partial["status"] == "harvested_partial_restart"
    assert partial["cp2k_program_ended"] is False
    assert partial["relaxation_completed_marker"] is False
    assert partial["optimization_steps"] == 41
    assert partial["relaxation_convergence_flags_all_yes"] is False

    completed_out = partial_out + "\n CELL OPTIMIZATION COMPLETED\n PROGRAM ENDED AT\n"
    completed = parse_relaxation_completion(completed_out, "CELL_OPT")
    assert completed["status"] == "cell_opt_converged"


def test_harvest_marks_restart_without_completion_as_partial(tmp_path):
    from pymatgen.core import Lattice, Structure
    from dft_validation.run_relaxation_validation import harvest

    results_dir = tmp_path / "relax"
    run_dir = results_dir / "ffc" / "runs" / "case_partial"
    run_dir.mkdir(parents=True)
    clean_cif = results_dir / "initial.cif"
    initial = Structure(
        Lattice.cubic(5.0),
        ["H", "H"],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        coords_are_cartesian=True,
    )
    initial.to(filename=str(clean_cif))
    (run_dir / "cp2k_run.out").write_text(
        """
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:            -1.000000000000
 OPTIMIZATION STEP:     3
 OPT| Maximum step size is converged                                          NO
"""
    )
    (run_dir / "cp2k_run-1_2.restart").write_text(
        """
&FORCE_EVAL
  &SUBSYS
    &CELL
      A 4.900000 0.000000 0.000000
      B 0.000000 4.900000 0.000000
      C 0.000000 0.000000 4.900000
    &END CELL
    &COORD
      SCALED .FALSE.
      H 0.000000 0.000000 0.000000
      H 0.000000 0.000000 0.740000
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
"""
    )
    records = [{
        "id": "case_partial",
        "topology": "ffc",
        "clean_cif": str(clean_cif),
        "n_atoms": 2,
        "status": "selected",
        "relax_run_type": "CELL_OPT",
    }]

    harvested = harvest(records, results_dir)
    assert harvested[0]["status"] == "harvested_partial_restart"
    assert harvested[0]["relaxed_extract_status"] == "cp2k_run-1_2.restart"
    assert harvested[0]["relaxed_cif"].endswith("case_partial.cif")


def test_cp2k_input_uses_numeric_latest_restart(tmp_path):
    import os
    from pymatgen.core import Lattice, Structure
    from dft_validation.cp2k_runner import generate_cp2k_input

    clean_cif = tmp_path / "initial.cif"
    Structure(
        Lattice.cubic(5.0),
        ["H", "H"],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        coords_are_cartesian=True,
    ).to(filename=str(clean_cif))
    for name in ["cp2k_run-1_9.restart", "cp2k_run-1_40.restart", "cp2k_run-1.restart"]:
        path = tmp_path / name
        path.write_text("&GLOBAL\n&END GLOBAL\n")
        os.utime(path, (1000, 1000))
    rec = {
        "id": "case_restart",
        "topology": "ffc",
        "clean_cif": str(clean_cif),
        "n_atoms": 2,
        "status": "selected",
    }

    inp = generate_cp2k_input(rec, str(tmp_path), "ffc", relax_mode_override="CELL_OPT")
    text = open(inp).read()
    assert "cp2k_run-1_40.restart" in text
    assert "cp2k_run-1_9.restart" not in text
