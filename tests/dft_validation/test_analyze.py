import sys
sys.path.insert(0, '.')
import pytest


def test_spearman_negative_trend():
    """Perfect inverse: higher rho → lower energy → rs < 0."""
    from dft_validation.analyze import _spearman_for_comp
    records = [
        {"rho": 1.0,   "dft_e_per_atom": -10.0, "sym_class": "High"},
        {"rho": 0.5,   "dft_e_per_atom": -9.5,  "sym_class": "High"},
        {"rho": 0.25,  "dft_e_per_atom": -9.0,  "sym_class": "Med"},
        {"rho": 0.125, "dft_e_per_atom": -8.5,  "sym_class": "Low"},
    ]
    rs, p, lowest_at_max_rho, n = _spearman_for_comp(records)
    assert rs < 0, f"Expected rs<0, got {rs}"
    assert lowest_at_max_rho is True
    assert n == 4


def test_retention_gate_passes_for_medium_high_minima():
    from dft_validation.analyze import get_fjh_gate
    comp_results = {
        "comp_a": {"n": 4, "lowest_rho": 0.5, "lowest_band": "Medium"},
        "comp_b": {"n": 4, "lowest_rho": 1.0, "lowest_band": "High"},
    }
    passed, reason = get_fjh_gate(comp_results)
    assert passed is True


def test_retention_gate_fails_on_low_minimum():
    from dft_validation.analyze import get_fjh_gate
    comp_results = {
        "comp_a": {"n": 4, "lowest_rho": 0.25, "lowest_band": "Low"},
        "comp_b": {"n": 4, "lowest_rho": 1.0, "lowest_band": "High"},
    }
    passed, reason = get_fjh_gate(comp_results)
    assert passed is False
    assert "comp_a" in reason


def test_retention_gate_pending_on_insufficient_converged_points():
    from dft_validation.analyze import get_fjh_gate
    comp_results = {
        "comp_a": {"n": 1, "lowest_rho": 0.5, "lowest_band": "Medium"},
        "comp_b": {"n": 4, "lowest_rho": 1.0, "lowest_band": "High"},
    }
    passed, reason = get_fjh_gate(comp_results)
    assert passed is None
    assert "insufficient" in reason
