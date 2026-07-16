import os
import sys
import random
import time
import csv
import numpy as np

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(current_dir)

try:
    from main import load_graph, calculate_topology_requirements
    from symmetry import SymmetryEngine
    from enumerator import MOFEnumerator
except ImportError:
    pass

def get_data_path(topology_name):
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")

def run_sampling_experiment_improved(topology_name, user_selected_ligands, num_samples=5000):
    print(f"\n>>> Running Baseline Experiment on [{topology_name}]")
    print(f"    Ligands: {user_selected_ligands}")

    graph_path = get_data_path(topology_name)
    if not os.path.exists(graph_path):
        print(f"Error: Graph file not found at {graph_path}")
        return [], [], []

    graph_data = load_graph(graph_path)
    all_ligands = sorted(set(lig for ligs in user_selected_ligands.values() for lig in ligs))
    node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
    topology_requirements = calculate_topology_requirements(graph_data)

    print("1. Computing Ground Truth (Burnside)...")
    sym = SymmetryEngine(graph_data, verbose=False)
    perms = sym.compute_automorphism_group()

    enum = MOFEnumerator(perms, graph_data['n_nodes'], verbose=False,
                         topology_requirements=topology_requirements,
                         user_selected_ligands=user_selected_ligands,
                         node_degrees=node_degrees)
    all_unique_structures = enum.enumerate(all_ligands)

    unique_orbits_map = {}
    high_sym_canonicals = set()

    for res in all_unique_structures:
        canon = tuple(res['coloring'])
        unique_orbits_map[canon] = res
        if res['stabilizer_size'] > 1:
            high_sym_canonicals.add(canon)

    total_unique_count = len(all_unique_structures)
    total_high_sym_targets = len(high_sym_canonicals)

    print(f"   -> Ground Truth: {total_unique_count} unique orbits.")
    print(f"   -> High-Sym Targets (Stab > 1): {total_high_sym_targets}")

    print(f"2. Random Sampling {num_samples} structures...")

    found_canonicals = set()
    found_history = []
    redundancy_history = []

    high_sym_coverage_history = []

    nodes = graph_data['n_nodes']
    atoms_pool = all_ligands

    start_time = time.time()

    for k in range(1, num_samples + 1):

        random_coloring = tuple(random.choices(atoms_pool, k=nodes))

        min_coloring = random_coloring
        for p in perms:
            mapped = tuple(random_coloring[i] for i in p)
            if mapped < min_coloring:
                min_coloring = mapped

        found_canonicals.add(min_coloring)

        unique_so_far = len(found_canonicals)
        found_history.append(unique_so_far)

        current_redundancy = 1.0 - (unique_so_far / k)
        redundancy_history.append(current_redundancy)

        current_high_sym_count = sum(1 for c in found_canonicals if c in high_sym_canonicals)

        if total_high_sym_targets > 0:
            current_coverage_pct = (current_high_sym_count / total_high_sym_targets) * 100.0
        else:
            current_coverage_pct = 0.0

        high_sym_coverage_history.append(current_coverage_pct)

    elapsed = time.time() - start_time

    final_high_sym_pct = high_sym_coverage_history[-1]

    global_cov = len(found_canonicals) / total_unique_count if total_unique_count > 0 else 0
    final_red = redundancy_history[-1]

    print(f"   -> Time: {elapsed:.2f}s")
    print(f"   -> Metric 1 (Global Coverage):   {len(found_canonicals)}/{total_unique_count} ({global_cov * 100:.2f}%)")
    print(f"   -> Metric 2 (High-Sym Coverage): {final_high_sym_pct:.2f}%")
    print(f"   -> Metric 3 (Redundancy):        {final_red * 100:.2f}%")

    return found_history, redundancy_history, high_sym_coverage_history

if __name__ == "__main__":

    SAMPLES = 10000
    test_cases = [
        # 4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"]}),

        # # 3c×1, 4c×2 → naive=1^16*2^12=4096
        # # 3c×2, 4c×1 → naive=2^16*1^12=65536
        # # 3c×1, 4c×3 → naive=1^16*3^12=531441
        ("fjh", {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"]}),

    ]
    TOPOLOGY, ATOMS = test_cases[0]

    atom_count = sum(len(v) for v in ATOMS.values())

    history_unique, history_red, history_high_cov = run_sampling_experiment_improved(TOPOLOGY, ATOMS,
                                                                                     num_samples=SAMPLES)

    if history_unique:

        output_dir = os.path.join(project_root, "results", "baseline")
        os.makedirs(output_dir, exist_ok=True)

        filename_base = f"experiment_data_{TOPOLOGY}_{SAMPLES}_{atom_count}"
        csv_path = os.path.join(output_dir, f"{filename_base}.csv")

        print(f"\n>>> Saving data to {csv_path}...")

        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)

            writer.writerow(["Sample_Index", "Unique_Found", "Redundancy_Rate", "High_Sym_Coverage_Pct"])
            for i in range(SAMPLES):
                writer.writerow([
                    i + 1,
                    history_unique[i],
                    history_red[i],
                    history_high_cov[i]
                ])

        print(">>> Plotting...")
        plt.figure(figsize=(4, 4))

        plt.subplot(1, 1, 1)
        plt.plot(history_high_cov, color='g', label='Stochastic Search')
        plt.xlabel('Number of Random Samples')

        plt.ylabel('Symmetric Target Coverage (%)')
        plt.title(f'{TOPOLOGY}, {atom_count} atoms')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        img_path = os.path.join(output_dir, f"{filename_base}.png")
        plt.savefig(img_path, dpi=300)
        print(f"   -> Plot Saved to {img_path}")