import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(current_dir)

from main import load_graph
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator
from verification import BridgeVerifier

def get_data_path(topology_name):
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")

def run_case(topology_name, atom_choices):
    print(f"[-] Processing {topology_name} with atoms: {atom_choices}...")

    graph_path = get_data_path(topology_name)

    if not os.path.exists(graph_path):
        print(f"Error: File not found at {graph_path}")
        return None

    graph_data = load_graph(graph_path)

    sym = SymmetryEngine(graph_data, verbose=False)
    perms = sym.compute_automorphism_group()

    enum = MOFEnumerator(perms, graph_data['n_nodes'], verbose=False)
    results = enum.enumerate(atom_choices)

    verifier = BridgeVerifier(len(perms))
    passed, stats = verifier.verify(results, graph_data['n_nodes'], atom_choices)
    return stats

if __name__ == "__main__":

    test_cases = [
        ("sod", ["6c_CdZn_1_Ch", "8c_Co_1_Ch"]),  # Case 1
        ("rht", ["6c_CdZn_1_Ch", "8c_Co_1_Ch"]),  # Case 2

        ("sod", ["6c_CdZn_1_Ch", "8c_Co_1_Ch", "6c_Al_1_Ch"]),
    ]

    print("=" * 80)
    print(f"{'Case':<25} | {'LHS (Constructive)':<20} | {'RHS (Theory)':<20} | {'Diff':<10}")
    print("-" * 80)

    for topo, atoms in test_cases:
        stats = run_case(topo, atoms)
        if stats:
            case_name = f"{topo}-{len(atoms)}types"
            print(f"{case_name:<25} | {stats['LHS']:<20.5f} | {stats['RHS']:<20.5f} | {stats['diff']:.2e}")

    print("=" * 80)