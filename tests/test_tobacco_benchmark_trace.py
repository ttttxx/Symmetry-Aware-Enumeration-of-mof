from pathlib import Path

from tests.benchmark_tobacco_reconstruction import (
    _degree_compatibility,
    _ligand_degree_map,
    _load_topology_context,
    _parse_case,
    _run_one_case,
)


def test_cds_4c_is_classified_before_assembly(tmp_path):
    case, error = _parse_case("cds_v1-4c_Cu_1_Ch_1B_4H_Ch.cif")

    assert error is None
    row = _run_one_case(case, tmp_path, max_cifs=3, max_results=5000)

    assert row["failure_class"] == "input_degree_incompatible"
    assert row["missing_topology_degrees"] == "3"
    assert row["n_generated_cifs"] == "0"
    assert not (tmp_path / "logs" / f"{Path(case.ref_file).stem}.log").exists()


def test_ffc_3c_4c_inputs_cover_topology_degrees():
    case, error = _parse_case(
        "ffc_v1-4c_Cu_1_Ch_v2-4c_Cu_1_Ch_v3-3c_B_Ch_1B_4H_Ch.cif"
    )

    assert error is None
    _, requirements, _ = _load_topology_context(case.topology)
    ligand_map = _ligand_degree_map(case.nodes)
    missing, extra = _degree_compatibility(requirements, ligand_map)

    assert requirements == {4: 3, 3: 4}
    assert sorted(ligand_map) == [3, 4]
    assert missing == []
    assert extra == []


def test_hwx_mixed_inputs_report_missing_4c_degree(tmp_path):
    case, error = _parse_case("hwx_v1-6c_Zn_1_Ch_v2-3c_B_Ch_1B_4H_Ch.cif")

    assert error is None
    row = _run_one_case(case, tmp_path, max_cifs=3, max_results=5000)

    assert row["failure_class"] == "input_degree_incompatible"
    assert row["missing_topology_degrees"] == "4"
    assert row["extra_ligand_degrees"] == "3"
    assert "degree_incompatible" in row["error_message"]
