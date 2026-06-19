import sys
sys.path.insert(0, '.')
import pytest
from pymatgen.core import Structure, Lattice
import numpy as np


def _make_struct(species, a=10.0):
    """Helper: cubic cell, atoms at fractional positions."""
    n = len(species)
    coords = [[i/n, 0.5, 0.5] for i in range(n)]
    return Structure(Lattice.cubic(a), species, coords)


def test_element_coverage_all_known():
    from dft_validation.sanitize import check_element_coverage
    struct = _make_struct(["C", "H", "O", "N", "Cd"])
    missing = check_element_coverage(struct)
    assert missing == [], f"Expected no missing, got {missing}"


def test_element_coverage_unknown_element():
    from dft_validation.sanitize import check_element_coverage
    struct = _make_struct(["C", "Xe"])  # Xe not in ELEMENT_KIND
    missing = check_element_coverage(struct)
    assert "Xe" in missing


def test_electron_parity_even():
    from dft_validation.sanitize import compute_electron_parity
    # C: 4, H: 1 → 4+1+1+4 = 10, even → multiplicity 1
    struct = _make_struct(["C", "H", "H", "C"])
    n_e, mult = compute_electron_parity(struct)
    assert n_e == 10
    assert mult == 1


def test_electron_parity_odd():
    from dft_validation.sanitize import compute_electron_parity
    # C:4, H:1 → 5, odd → multiplicity 2
    struct = _make_struct(["C", "H"])
    n_e, mult = compute_electron_parity(struct)
    assert n_e == 5
    assert mult == 2


def test_min_distance_clean_struct():
    """A well-spaced structure reports min_dist above floor."""
    from dft_validation.sanitize import min_nonbond_distance
    struct = _make_struct(["C", "C", "C"], a=15.0)
    d = min_nonbond_distance(struct)
    assert d > 2.0


def test_xtb_runs_on_cd_structure():
    """xTB FIRE runs without crash on a real assembled Cd CIF (integration)."""
    import glob
    cifs = glob.glob("results/dft_validation/ffc/init/*comp_a*.cif")
    if not cifs:
        pytest.skip("Phase 1 output not present — run --phase 1 first")
    from dft_validation.sanitize import preclean_xtb
    from pymatgen.core import Structure
    struct = Structure.from_file(cifs[0])
    result = preclean_xtb(struct)
    assert result is not None, "xTB returned None"
    assert hasattr(result, "sites"), "Result is not a Structure"
