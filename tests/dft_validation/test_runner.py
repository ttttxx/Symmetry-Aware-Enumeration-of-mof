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
    rec = next(r for r in recs if r.get("status") in ("clean", "preclean_done")
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
    rec = next(r for r in recs if r.get("status") in ("clean", "preclean_done"))
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
