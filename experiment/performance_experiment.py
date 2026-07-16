import os
import sys
import time
import csv
import math
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(current_dir)
sys.path.append(project_root)

from main import load_graph, calculate_topology_requirements
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator

def get_data_path(topology_name):
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")

def benchmark_case(topology_name, user_selected_ligands):
    path = get_data_path(topology_name)
    if not os.path.exists(path):
        print(f"Warning: {path} not found.")
        return None

    graph_data = load_graph(path)
    n_nodes = graph_data['n_nodes']
    node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
    topology_requirements = calculate_topology_requirements(graph_data)

    all_ligands = sorted(set(
        lig for ligs in user_selected_ligands.values() for lig in ligs
    ))
    total_atom_types = len(all_ligands)

    label_parts = []
    for deg in sorted(user_selected_ligands.keys()):
        label_parts.append(f"{deg}c×{len(user_selected_ligands[deg])}")
    ligand_label = ",".join(label_parts)

    naive_space = 1
    for deg, count in topology_requirements.items():
        n_choices = len(user_selected_ligands.get(deg, []))
        if n_choices > 0:
            naive_space *= (n_choices ** count)
    naive_space_log10 = math.log10(naive_space) if naive_space > 0 else 0.0

    print(f"   Benchmarking {topology_name} | ligands: {ligand_label} "
          f"| nodes: {n_nodes} | deg-dist: {topology_requirements} | "
          f"naive space: 10^{naive_space_log10:.1f}")

    t0 = time.time()
    sym = SymmetryEngine(graph_data, verbose=False)
    perms = sym.compute_automorphism_group()
    t_sym = time.time() - t0
    group_size = len(perms)

    t1 = time.time()
    enumerator = MOFEnumerator(
        perms, n_nodes, verbose=False,
        topology_requirements=topology_requirements,
        user_selected_ligands=user_selected_ligands,
        node_degrees=node_degrees
    )
    results = enumerator.enumerate(all_ligands)
    t_enum = time.time() - t1

    n_unique = len(results)

    return {
        "Topology": topology_name,
        "Ligand_Label": ligand_label,
        "Total_Atom_Types": total_atom_types,
        "Nodes": n_nodes,
        "Degree_Dist": str(topology_requirements),
        "Group_Size": group_size,
        "Naive_Log10": round(naive_space_log10, 2),
        "Unique_Found": n_unique,
        "Time_Sym": round(t_sym, 4),
        "Time_Enum": round(t_enum, 4),
        "Time_Total": round(t_sym + t_enum, 4),
    }

def run_performance_experiment():
    print(">>> Starting Complexity & Performance Benchmark (Degree-Aware)...")

    # Naive Space = Π |ligands[deg]|^{count[deg]}

    test_cases = [

        # 4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"]}),

        # # 3c×1, 4c×2 → naive=1^16*2^12=4096
        # # 3c×2, 4c×1 → naive=2^16*1^12=65536
        # # 3c×1, 4c×3 → naive=1^16*3^12=531441

        ("vba", {4: ["4c_Cu_1_Ch", "4c_B_Ch"]}),
        ("vba", {4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"]}),
    ]

    results = []

    output_dir = os.path.join(project_root, "results", "performance")
    os.makedirs(output_dir, exist_ok=True)

    for topo, ligands in test_cases:
        data = benchmark_case(topo, ligands)
        if data:
            results.append(data)

    csv_path = os.path.join(output_dir, "performance_data.csv")
    fieldnames = [
        "Topology", "Ligand_Label", "Total_Atom_Types", "Nodes",
        "Degree_Dist", "Group_Size", "Naive_Log10",
        "Unique_Found", "Time_Sym", "Time_Enum", "Time_Total"
    ]
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f">>> Data saved to {csv_path}")
    plot_performance(results, output_dir)

def plot_performance(data, output_dir):
    print(">>> Plotting Performance Chart...")

    labels = [f"{d['Topology']}\n{d['Ligand_Label']}" for d in data]
    times = [d['Time_Total'] for d in data]
    naive_logs = [d['Naive_Log10'] for d in data]
    unique_counts = [d['Unique_Found'] for d in data]

    fig, ax1 = plt.subplots(figsize=(14, 7))

    color_bar = 'tab:blue'
    ax1.set_xlabel('Test Case (Topology — Ligand Configuration)', fontweight='bold', fontsize=11)
    ax1.set_ylabel('Actual Runtime (seconds)', color=color_bar, fontweight='bold')
    bars = ax1.bar(labels, times, color=color_bar, alpha=0.6, width=0.5, label='Actual Time')
    ax1.tick_params(axis='y', labelcolor=color_bar)
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, fontsize=9)

    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2., height,
                 f'{height:.2f}s\n({unique_counts[i]} uniq)',
                 ha='center', va='bottom', color='black', fontsize=8)

    ax2 = ax1.twinx()
    color_line = 'tab:red'
    ax2.set_ylabel('Naive Search Space ($10^N$)', color=color_line, fontweight='bold')
    ax2.plot(range(len(labels)), naive_logs, color=color_line, marker='o',
             linewidth=2, linestyle='--', markersize=8, label='Naive Space')
    ax2.tick_params(axis='y', labelcolor=color_line)

    for i, val in enumerate(naive_logs):
        ax2.text(i, val + 0.3, f'$10^{{{val:.1f}}}$',
                 ha='center', va='bottom', color='red', fontsize=8)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)

    plt.title('Performance vs. Combinatorial Explosion (Degree-Aware Enumeration)',
              fontsize=14, fontweight='bold')
    fig.tight_layout()

    plot_path = os.path.join(output_dir, "complexity_performance.png")
    plt.savefig(plot_path, dpi=300)
    print(f"   -> Plot saved to {plot_path}")

if __name__ == "__main__":
    run_performance_experiment()
