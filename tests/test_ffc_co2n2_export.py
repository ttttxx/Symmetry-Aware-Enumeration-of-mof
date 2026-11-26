import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, ".")


def test_node_screening_suffixes_ethene_as_2c():
    from adsorption_screening.export_ffc_co2n2_candidates import screen_nodes

    infos = {info.name: info for info in screen_nodes()}

    assert infos["4c_ethene_Ch"].inferred_connectivity == 2
    assert infos["4c_ethene_Ch"].snapshot_name == "4c_ethene_Ch_2c.cif"
    assert not infos["4c_ethene_Ch"].selected_for_batch


def test_selected_node_pool_excludes_risky_or_mismatched_nodes():
    from adsorption_screening.export_ffc_co2n2_candidates import screen_nodes, selected_atom_choices

    infos = screen_nodes()
    choices = selected_atom_choices(infos)
    selected = {name for names in choices.values() for name in names}

    assert "3c_uranyl" not in selected
    assert "3c_uranyl_Ch" not in selected
    assert "4c_1anC_Ch" not in selected
    assert "4c_ethene_Ch" not in selected
    assert "4c_Zn_tet" not in selected
    assert "4c_pyrene_NU_Ch" not in selected
    assert choices[3] == ["3c_B_Ch", "3c_Cu_imadazoline_Ch", "3c_pyrrole_Ch", "3c_triazine_Ch"]
    assert len(choices[4]) == 10


def test_dry_run_selection_returns_exact_batch_counts():
    from adsorption_screening.export_ffc_co2n2_candidates import export_batch

    rows = export_batch(dry_run=True)
    buckets = Counter(row["selection_bucket"] for row in rows)

    assert len(rows) == 200
    assert buckets == {
        "symmetry_priority": 140,
        "descriptor_diversity": 40,
        "low_sym_control": 20,
    }


def test_dry_run_manifest_minimum_fields_are_present():
    from adsorption_screening.export_ffc_co2n2_candidates import export_batch

    rows = export_batch(dry_run=True)
    core_required = {"candidate_id", "topology", "composition", "rho", "sym_class", "cif_path"}
    validation_required = {
        "validation_passed",
        "clash_passed",
        "topology_passed",
        "n_clashes",
        "validation_json",
    }
    required = core_required | validation_required

    for row in rows:
        assert required <= set(row)
        for key in core_required:
            assert row[key] not in (None, "")
        assert row["topology"] == "ffc"


def test_dry_run_manifest_validation_fields_are_present():
    from adsorption_screening.export_ffc_co2n2_candidates import export_batch

    rows = export_batch(dry_run=True)
    required = {
        "candidate_id",
        "validation_passed",
        "clash_passed",
        "topology_passed",
        "n_clashes",
        "validation_json",
    }

    for row in rows:
        assert required <= set(row)


def test_validation_manifest_fields_are_json_serializable(tmp_path):
    from pymatgen.core import Lattice, Structure

    from adsorption_screening.export_ffc_co2n2_candidates import _manifest_row

    record = {
        "composition": "3c_B_Ch:1",
        "rho": 0.5,
        "sym_class": "High",
        "selection_bucket": "symmetry_priority",
        "source_enum_id": 1,
        "stabilizer_size": 4,
        "gamma_size": 8,
        "coloring": ["3c_B_Ch"],
    }
    node = type(
        "NodeInfoLike",
        (),
        {"path": "nodes/3c_B_Ch.cif", "family_tags": ["hydrocarbon"]},
    )()
    structure = Structure(Lattice.cubic(10), ["C"], [[0, 0, 0]])
    metadata = {
        "validation_passed": False,
        "detected_sg": "P1",
        "validation": {
            "passed": False,
            "clash": {
                "passed": False,
                "n_clashes": 2,
                "worst_clash": {"ratio": 0.42},
                "message": "clash",
            },
            "topology": {"passed": True, "message": "ok"},
            "symmetry": {"passed": True, "message": "ok", "detected_sg": "P1"},
        },
    }

    row = _manifest_row(record, "demo", tmp_path / "demo.cif", structure, {"3c_B_Ch": node}, metadata)

    assert row["validation_passed"] is False
    assert row["clash_passed"] is False
    assert row["topology_passed"] is True
    assert row["n_clashes"] == 2
    assert row["worst_clash_ratio"] == 0.42
    assert row["validation_warning"] == "clash"
    assert json.loads(row["validation_json"])["clash"]["n_clashes"] == 2


def test_real_export_to_tmp_path_writes_reloadable_cifs(tmp_path):
    from pymatgen.core import Structure

    from adsorption_screening.export_ffc_co2n2_candidates import export_batch

    rows = export_batch(out_root=tmp_path)
    batch_dir = tmp_path / "batch_001"
    manifest_path = batch_dir / "manifest.csv"
    nodes_path = batch_dir / "nodes_screening.csv"
    snapshot_path = batch_dir / "nodes_snapshot" / "4c_ethene_Ch_2c.cif"

    assert len(rows) == 200
    assert manifest_path.exists()
    assert nodes_path.exists()
    assert snapshot_path.exists()

    with manifest_path.open(newline="", encoding="utf-8") as handle:
        manifest_rows = list(csv.DictReader(handle))

    assert len(manifest_rows) == 200
    assert {row["topology"] for row in manifest_rows} == {"ffc"}
    for row in manifest_rows:
        assert row["validation_passed"] in {"True", "False"}
        assert row["clash_passed"] in {"True", "False"}
        assert row["topology_passed"] in {"True", "False"}
        assert row["n_clashes"] != ""
        assert isinstance(json.loads(row["validation_json"]), dict)
        cif_path = Path(row["cif_path"])
        assert cif_path.exists()
        structure = Structure.from_file(cif_path)
        assert len(structure) == int(row["n_atoms"])
        assert len(structure) <= 450


def test_clean_pool_node_set_has_large_enumeration_space():
    from adsorption_screening.export_ffc_co2n2_clean_pool import CLEAN_BIASED_NODE_CHOICES
    from adsorption_screening.export_ffc_co2n2_candidates import _enumerate_ffc

    records, _graph_data, gamma = _enumerate_ffc(CLEAN_BIASED_NODE_CHOICES)
    classes = Counter(record["sym_class"] for record in records)

    assert gamma == 8
    assert len(records) >= 5000
    assert classes["High"] > 0
    assert classes["Med"] > 0
    assert classes["Low"] > 0


def test_clean_pool_export_writes_only_strict_clean_manifest_rows(tmp_path):
    from pymatgen.core import Structure

    from adsorption_screening.export_ffc_co2n2_clean_pool import export_clean_pool

    rows = export_clean_pool(
        batch_id="smoke_clean_pool",
        target_clean=3,
        out_root=tmp_path,
        progress_every=0,
    )
    batch_dir = tmp_path / "smoke_clean_pool"

    assert len(rows) == 3
    assert (batch_dir / "manifest.csv").exists()
    assert (batch_dir / "rejected_manifest.csv").exists()
    assert (batch_dir / "summary.json").exists()

    with (batch_dir / "manifest.csv").open(newline="", encoding="utf-8") as handle:
        manifest_rows = list(csv.DictReader(handle))

    assert len(manifest_rows) == 3
    for row in manifest_rows:
        assert row["validation_passed"] == "True"
        assert row["clash_passed"] == "True"
        assert row["topology_passed"] == "True"
        cif_path = Path(row["cif_path"])
        assert cif_path.exists()
        structure = Structure.from_file(cif_path)
        assert len(structure) == int(row["n_atoms"])
