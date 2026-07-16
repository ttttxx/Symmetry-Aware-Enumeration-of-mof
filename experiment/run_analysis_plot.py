import os
import sys
import json
import csv
import matplotlib.pyplot as plt
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(project_root)

try:
    from main import load_graph
except ImportError:

    def load_graph(path):
        with open(path, 'r') as f:
            return json.load(f)

results_root = os.path.join(project_root, "results", "structure_selection")

BOND_ENERGIES = {
    "same_type": 0.10,
    "diff_type": -0.20
}

def get_data_path(topology_name):

    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")

def calculate_ising_energy(structure_data, graph_data):
    coloring = structure_data['coloring']

    edges_raw = None
    adj_list = None

    if isinstance(graph_data, list):
        adj_list = graph_data
    elif 'edges' in graph_data:
        edges_raw = graph_data['edges']
    elif 'adjacency' in graph_data:
        adj_list = graph_data['adjacency']

    total_energy = 0.0
    valid_edges_count = 0

    if edges_raw is not None:
        for edge in edges_raw:
            u, v = None, None

            if isinstance(edge, dict):

                u = edge.get('source', edge.get('u'))
                v = edge.get('target', edge.get('v'))

                if u is None: u = edge.get('node1')
                if v is None: v = edge.get('node2')
            elif isinstance(edge, (list, tuple)):

                if len(edge) >= 2:
                    u = edge[0]
                    v = edge[1]

            if u is None or v is None:
                continue

            if u >= len(coloring) or v >= len(coloring):
                continue

            type_u = coloring[u]
            type_v = coloring[v]

            if type_u == type_v:
                total_energy += BOND_ENERGIES['same_type']
            else:
                total_energy += BOND_ENERGIES['diff_type']
            valid_edges_count += 1

        if valid_edges_count > 0:
            return round(total_energy, 4)

    if adj_list is not None:
        for u, neighbors in enumerate(adj_list):
            type_u = coloring[u]
            for val in neighbors:

                v = None
                if isinstance(val, (int, np.integer)):
                    v = int(val)
                elif isinstance(val, (list, tuple)):
                    v = int(val[0])
                elif isinstance(val, dict):
                    v = val.get('id', val.get('node', None))

                if v is None: continue

                if u < v:
                    type_v = coloring[v]
                    if type_u == type_v:
                        total_energy += BOND_ENERGIES['same_type']
                    else:
                        total_energy += BOND_ENERGIES['diff_type']
        return round(total_energy, 4)

    print("[Warning] Could not parse edge connectivity from graph data.")
    return 0.0

def run_analysis(topology_name):

    target_dir = os.path.join(results_root, f"{topology_name}_selection")
    topology_file = get_data_path(topology_name)

    if not os.path.exists(target_dir):
        print(f"[Error] Directory not found: {target_dir}")
        return

    print(f"[-] Loading Topology Graph from: {topology_file}")
    if not os.path.exists(topology_file):
        print(f"[Error] Topology file missing: {topology_file}")
        return

    graph_data = load_graph(topology_file)

    print(f"[-] Reading structures from: {target_dir}")
    json_files = [f for f in os.listdir(target_dir) if f.endswith(".json")]

    if not json_files:
        print("[Error] No JSON files found in selection folder.")
        return

    data_points = []

    print(f"[-] Calculating energy using Bond Valence Model (Ising)...")
    for filename in json_files:
        filepath = os.path.join(target_dir, filename)

        with open(filepath, 'r') as f:
            struct_data = json.load(f)

        energy = calculate_ising_energy(struct_data, graph_data)

        struct_data['properties']['energy_ev'] = energy

        label = struct_data['label']
        if "High" in label:
            stype = "High"
        elif "Med" in label:
            stype = "Med"
        else:
            stype = "Low"

        data_points.append({
            "filename": filename,
            "label": label,
            "sym_order": struct_data['symmetry_order'],
            "energy": energy,
            "type": stype
        })

    data_points.sort(key=lambda x: x['energy'])

    csv_path = os.path.join(results_root, f"{topology_name}_analysis_data.csv")
    save_csv(data_points, csv_path)

    plot_scatter(data_points, topology_name)
    print("Done.")

def save_csv(data, filepath):
    headers = ["Filename", "Label", "Symmetry_Order", "Energy_eV", "Type"]

    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in data:
            writer.writerow([
                row['filename'],
                row['label'],
                row['sym_order'],
                row['energy'],
                row['type']
            ])
    print(f"[Output] Data saved to CSV: {filepath}")

def plot_scatter(data, topology_name):
    plt.figure(figsize=(10, 7))

    style_map = {
        'High': {'color': '#d62728', 'marker': 'o', 'label': 'High Symmetry (Top-10)'},  # Red
        'Med': {'color': '#1f77b4', 'marker': '^', 'label': 'Medium Symmetry'},  # Blue
        'Low': {'color': '#7f7f7f', 'marker': 'x', 'label': 'Asymmetric (Low-10)'}  # Grey
    }

    groups = {'High': [], 'Med': [], 'Low': []}
    for p in data:
        groups[p['type']].append(p)

    for g_name in ['High', 'Med', 'Low']:
        points = groups[g_name]
        if not points: continue

        x = [p['sym_order'] for p in points]
        y = [p['energy'] for p in points]
        style = style_map[g_name]

        if style['marker'] == 'x':
            plt.scatter(x, y, c=style['color'], marker=style['marker'],
                        label=style['label'], s=120, alpha=0.8, linewidths=2)
        else:
            plt.scatter(x, y, c=style['color'], marker=style['marker'],
                        label=style['label'], s=120, alpha=0.8, edgecolors='k')

    plt.title(f'Symmetry vs Stability: {topology_name.upper()} Topology', fontsize=16, fontweight='bold')
    plt.xlabel('Symmetry Order (Size of Stabilizer Group)', fontsize=14)
    plt.ylabel('Configurational Energy (eV) [Ising Model]', fontsize=14)

    plt.xscale('log')

    plt.grid(True, linestyle='--', alpha=0.4, which="both")
    plt.legend(fontsize=12)
    plt.tight_layout()

    output_img = os.path.join(results_root, f"{topology_name}_analysis.png")
    plt.savefig(output_img, dpi=300)
    print(f"[Output] Plot saved to PNG: {output_img}")

if __name__ == "__main__":
    run_analysis("sod")