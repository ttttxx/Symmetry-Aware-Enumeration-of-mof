import sys, os
sys.path.insert(0, '.')
import pytest
from collections import Counter


def _apply_cd_subst(atom_choices):
    """Mirror of select_structures._apply_cd_substitution."""
    result = {}
    for deg, ligs in atom_choices.items():
        result[deg] = ["4c_Cd_1_Ch" if l == "4c_Cu_1_Ch" else l for l in ligs]
    return result


def test_comp_key_matches_expected():
    from dft_validation.select_structures import comp_key
    coloring = ["3c_B_Ch", "3c_triazine_Ch", "3c_B_Ch", "3c_triazine_Ch",
                "4c_Cd_1_Ch", "4c_B_Ch", "4c_B_Ch"]
    key = comp_key(coloring)
    assert key == {"3c_B_Ch": 2, "3c_triazine_Ch": 2, "4c_Cd_1_Ch": 1, "4c_B_Ch": 2}


def test_classify_stabs_ffc():
    from dft_validation.select_structures import classify_sym_classes
    # ffc stabs = {8,4,2,1} → high={8,4}, med={2}, low={1}
    fake = [{"stabilizer_size": s} for s in [8, 8, 4, 4, 2, 2, 1, 1]]
    classified = classify_sym_classes(fake)
    classes = [r["sym_class"] for r in classified]
    assert classes.count("High") == 4
    assert classes.count("Med")  == 2
    assert classes.count("Low")  == 2


def test_ffc_gives_10_structures():
    """Integration: real enumeration → exactly 10 ffc Cd structures (5 A + 5 B)."""
    from dft_validation.select_structures import enumerate_ffc
    from dft_validation import config
    target = next(t for t in config.VALIDATION_TARGETS if t["topology"] == "ffc")
    records = enumerate_ffc(target)
    assert len(records) == 10, f"Expected 10, got {len(records)}"
    labels = [r["composition_label"] for r in records]
    assert labels.count("comp_a") == 5
    assert labels.count("comp_b") == 5


def test_ffc_sym_classes_present():
    from dft_validation.select_structures import enumerate_ffc
    from dft_validation import config
    target = next(t for t in config.VALIDATION_TARGETS if t["topology"] == "ffc")
    records = enumerate_ffc(target)
    for label in ["comp_a", "comp_b"]:
        subset = [r for r in records if r["composition_label"] == label]
        classes = {r["sym_class"] for r in subset}
        assert "High" in classes, f"{label} missing High"
        assert "Med"  in classes, f"{label} missing Med"
        assert "Low"  in classes, f"{label} missing Low"
