import os
import sys
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from enumerator import IDENTITY_ORBIT_LIMIT, IDENTITY_NEW_LIMIT

def test_constants_exist_and_positive():
    assert IDENTITY_ORBIT_LIMIT > 0
    assert IDENTITY_NEW_LIMIT > 0
    assert IDENTITY_ORBIT_LIMIT > IDENTITY_NEW_LIMIT

from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator

def _make_ffc_enumerator():
    """返回 ffc 拓扑的枚举器（7节点，|G|=8，小而可控）"""
    base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    graph_data = load_graph(os.path.join(base, "data", "ffc", "ffc_graph.json"))
    graph_data = enrich_graph_data_with_topology(graph_data)
    node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
    topo_req = calculate_topology_requirements(graph_data)
    ligands = {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]}
    all_ligs = sorted(set(l for ls in ligands.values() for l in ls))
    perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
    return MOFEnumerator(perms, graph_data['n_nodes'], verbose=False,
                         topology_requirements=topo_req,
                         user_selected_ligands=ligands,
                         node_degrees=node_degrees), all_ligs

def test_budget_limits_new_results():
    """当 identity_budget=2 时，base case 不应新增超过 2 个结构"""
    enumerator, all_ligs = _make_ffc_enumerator()
    identity_orbits = [[i] for i in range(enumerator.n)]
    initial_active = list(range(1, len(enumerator.perms)))
    current_coloring = [None] * enumerator.n
    budget = 2
    start = len(enumerator.classes)
    enumerator._dfs_lattice_pruning(
        0, identity_orbits, all_ligs, current_coloring, initial_active,
        identity_start=start, identity_budget=budget
    )
    added = len(enumerator.classes) - start
    assert 1 <= added <= budget, f"Expected 1..{budget} new results, got {added}"


def test_small_input_full_enumeration():
    """ffc + {3c×2, 4c×2}：naive space 小，Identity 应完整枚举，Burnside SUCCESS"""
    enumerator, all_ligs = _make_ffc_enumerator()
    results = enumerator.enumerate(all_ligs)
    assert len(results) > 0
    lhs = sum(1.0 / r['stabilizer_size'] for r in results)
    naive = 1
    topo = enumerator.topology_requirements
    for deg, cnt in topo.items():
        n_ch = len(enumerator.user_selected_ligands.get(deg, all_ligs))
        naive *= n_ch ** cnt
    rhs = naive / len(enumerator.perms)
    assert abs(lhs - rhs) < 1e-6, f"Burnside mismatch: lhs={lhs:.6f} rhs={rhs:.6f}"

def test_large_input_partial_enumeration_via_low_limit():
    """调低阈值模拟大输入：Identity 应限额枚举"""
    import enumerator as enum_module
    original_limit = enum_module.IDENTITY_ORBIT_LIMIT
    original_new = enum_module.IDENTITY_NEW_LIMIT
    try:
        enum_module.IDENTITY_ORBIT_LIMIT = 1
        enum_module.IDENTITY_NEW_LIMIT = 3
        en, ligs = _make_ffc_enumerator()
        results = en.enumerate(ligs)
        # Partial mode should still return some results (not zero) and not crash
        assert len(results) > 0, "Partial Identity enumeration should still return some results"
        # With IDENTITY_NEW_LIMIT=3, Identity adds at most 3 new structures.
        # Higher-symmetry tasks may have already found some, so total can be more than 3.
        # But the run must complete without hanging.
    finally:
        enum_module.IDENTITY_ORBIT_LIMIT = original_limit
        enum_module.IDENTITY_NEW_LIMIT = original_new

def test_identity_skipped_prints_warning(capsys):
    """调低阈值后，verbose 模式应打印限额警告"""
    import enumerator as enum_module
    original_limit = enum_module.IDENTITY_ORBIT_LIMIT
    try:
        enum_module.IDENTITY_ORBIT_LIMIT = 1
        base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        graph_data = load_graph(os.path.join(base, "data", "ffc", "ffc_graph.json"))
        graph_data = enrich_graph_data_with_topology(graph_data)
        node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
        topo_req = calculate_topology_requirements(graph_data)
        ligands = {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]}
        all_ligs = sorted(set(l for ls in ligands.values() for l in ls))
        perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
        en = MOFEnumerator(perms, graph_data['n_nodes'], verbose=True,
                           topology_requirements=topo_req,
                           user_selected_ligands=ligands,
                           node_degrees=node_degrees)
        en.enumerate(all_ligs)
        out = capsys.readouterr().out
        assert "Partial" in out or "partial" in out or "budget" in out or "Budget" in out
    finally:
        enum_module.IDENTITY_ORBIT_LIMIT = original_limit
