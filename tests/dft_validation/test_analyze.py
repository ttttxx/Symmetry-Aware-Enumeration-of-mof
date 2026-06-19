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
    rs, p, lowest_in_high = _spearman_for_comp(records)
    assert rs < 0, f"Expected rs<0, got {rs}"
    assert lowest_in_high is True


def test_fjh_gate_passes():
    from dft_validation.analyze import get_fjh_gate
    comp_results = {
        "comp_a": {"spearman_rs": -0.9, "p_val": 0.05, "lowest_in_high": True},
        "comp_b": {"spearman_rs": -0.8, "p_val": 0.08, "lowest_in_high": True},
    }
    passed, reason = get_fjh_gate(comp_results)
    assert passed is True


def test_fjh_gate_fails_on_not_in_high():
    from dft_validation.analyze import get_fjh_gate
    comp_results = {
        "comp_a": {"spearman_rs": -0.9, "p_val": 0.05, "lowest_in_high": False},
        "comp_b": {"spearman_rs": -0.8, "p_val": 0.08, "lowest_in_high": True},
    }
    passed, reason = get_fjh_gate(comp_results)
    assert passed is False
    assert "comp_a" in reason


def test_fjh_gate_fails_on_positive_rs():
    from dft_validation.analyze import get_fjh_gate
    comp_results = {
        "comp_a": {"spearman_rs": 0.3,  "p_val": 0.5,  "lowest_in_high": True},
        "comp_b": {"spearman_rs": -0.8, "p_val": 0.08, "lowest_in_high": True},
    }
    passed, reason = get_fjh_gate(comp_results)
    assert passed is False
