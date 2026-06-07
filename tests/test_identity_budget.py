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
    assert added <= budget, f"Expected ≤{budget} new results, got {added}"
