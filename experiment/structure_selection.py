import os
import sys
import json
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))

project_root = os.path.dirname(current_dir)

sys.path.append(project_root)

from main import load_graph
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator

def get_data_path(topology_name):
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")

def ensure_output_dir(subfolder):

    out_path = os.path.join(project_root, "results", "structure_selection", subfolder)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    return out_path

def select_and_save_structures(topology_name, atom_choices):
    print(f"\n>>> Running Structure Selection for Analysis: {topology_name} {atom_choices}")

    graph_path = get_data_path(topology_name)
    if not os.path.exists(graph_path):
        print(f"Error: Data file not found at {graph_path}")
        return

    graph_data = load_graph(graph_path)
    sym = SymmetryEngine(graph_data, verbose=False)
    perms = sym.compute_automorphism_group()

    enum = MOFEnumerator(perms, graph_data['n_nodes'], verbose=False)
    all_results = enum.enumerate(atom_choices)

    print(f"   -> Total unique structures found: {len(all_results)}")

    sorted_structs = sorted(all_results, key=lambda x: x['stabilizer_size'], reverse=True)

    selected_samples = []

    top_5 = sorted_structs[:10]
    for i, res in enumerate(top_5):
        res['_label'] = f"HighSym_{i + 1}"
        selected_samples.append(res)

    low_5 = sorted_structs[-10:]

    for i, res in enumerate(low_5):
        res['_label'] = f"LowSym_{i + 1}"
        selected_samples.append(res)

    mid_idx = len(sorted_structs) // 2
    mid_start = max(0, mid_idx - 2)
    mid_10 = sorted_structs[mid_start: mid_start + 10]
    for i, res in enumerate(mid_10):
        res['_label'] = f"MedSym_{i + 1}"
        selected_samples.append(res)

    output_dir = ensure_output_dir(f"{topology_name}_selection")

    print(f"   -> Saving {len(selected_samples)} selected structures to {output_dir}...")

    summary_list = []

    for res in selected_samples:
        label = res['_label']

        filename = f"{topology_name}_{label}.json"
        save_path = os.path.join(output_dir, filename)

        save_data = {
            "topology": topology_name,
            "label": label,
            "symmetry_order": res['stabilizer_size'],
            "coloring": res['coloring'],
            "atom_choices": atom_choices,

            "properties": {
                "energy_ev": None,
                "stability_score": None
            }
        }

        with open(save_path, 'w') as f:
            json.dump(save_data, f, indent=4)

        summary_list.append({
            "Label": label,
            "Sym_Order": res['stabilizer_size'],
            "Coloring_Snippet": str(res['coloring'][:6]) + "..."
        })

    print(f"\n{'-' * 60}")
    print(f"{'Label':<15} | {'Sym Order':<10} | {'Structure Snippet'}")
    print(f"{'-' * 60}")
    for item in summary_list:
        print(f"{item['Label']:<15} | {item['Sym_Order']:<10} | {item['Coloring_Snippet']}")
    print(f"{'-' * 60}")
    print(f"Done. You can now use files in {output_dir} for stability analysis.\n")

if __name__ == "__main__":

    select_and_save_structures("sod", ["6c_CdZn_1_Ch", "8c_Co_1_Ch"])

