import os
import sys
import csv
import random
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from chgnet.model import CHGNet
from pymatgen.core import Structure

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(current_dir)
sys.path.append(project_root)

from main import load_graph
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator

from assembler import StructureAssembler

print(">>> Loading CHGNet Universal Potential (Replacing Ising Model)...")
chgnet_model = CHGNet.load()

STYLE_MAP = {
    'High': {'color': '#d62728', 'marker': 'o'},
    'Med': {'color': '#1f77b4', 'marker': '^'},
    'Low': {'color': '#95a5a6', 'marker': 'D'}
}

def get_data_path(topology_name):
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")

def evaluate_physical_energy(assembler, coloring, graph_data, temp_cif_path="temp.cif", max_atoms=2500):
    try:
        structure_map = {str(i): color for i, color in enumerate(coloring)}
        assembler.assemble_and_save(
            filename=temp_cif_path,
            structure_map=structure_map,
            graph_data=graph_data
        )
        struct = Structure.from_file(temp_cif_path)

        atom_count = len(struct)
        if atom_count > max_atoms:
            print(f"      [Skip] 结构过大 (含 {atom_count} 个原子 > 限制 {max_atoms}), 自动跳过以防内存崩溃。")
            return None

        prediction = chgnet_model.predict_structure(struct)
        energy_per_atom = prediction['e']
        return energy_per_atom
    except Exception as e:
        print(f"      [Warning] Reconstruction or Evaluation failed: {e}")
        return None

def run_physical_experiment(experiments, samples_per_category=20):
    print(f">>> Starting Thermodynamic Stability Validation (CHGNet Prediction)")

    output_dir = os.path.join(project_root, "results", "physical_validation")
    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, f"real_physics_stability_data_{samples_per_category}.csv")

    os.chdir(project_root)
    assembler = StructureAssembler(output_dir=output_dir)
    all_datasets = []

    with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Topology", "Atom_Count", "Label", "Symmetry_Index_Rho", "Energy_eV_per_atom", "Category"])

        for topology_name, atom_choices in experiments:
            label = f"{topology_name}_{len(atom_choices)}"
            print(f"\n--- Processing: {label} ---")

            graph_path = get_data_path(topology_name)
            if not os.path.exists(graph_path):
                print(f"Skipping {label}: File not found.")
                continue

            graph_data = load_graph(graph_path)

            sym = SymmetryEngine(graph_data, verbose=False)
            perms = sym.compute_automorphism_group()

            gamma_size = len(perms)
            print(f"   Automorphism Group Size |\Gamma|: {gamma_size}")

            enum = MOFEnumerator(perms, graph_data['n_nodes'], verbose=False)
            results = enum.enumerate(atom_choices)
            print(f"   Total unique structures enumerated: {len(results)}")

            all_stabs = sorted(list(set(r['stabilizer_size'] for r in results)), reverse=True)
            high_vals, med_vals = set(), set()

            if len(all_stabs) <= 2:
                if all_stabs[0] > 1: high_vals.add(all_stabs[0])
            elif len(all_stabs) == 3:
                high_vals.add(all_stabs[0])
                if all_stabs[1] > 1: med_vals.add(all_stabs[1])
            else:
                high_vals.update([all_stabs[0], all_stabs[1]])
                for v in all_stabs[2:]:
                    if v > 1: med_vals.add(v)

            high_sym = [r for r in results if r['stabilizer_size'] in high_vals]
            med_sym = [r for r in results if r['stabilizer_size'] in med_vals]
            low_sym = [r for r in results if r['stabilizer_size'] == 1]

            max_high_eval = 300
            high_sym_sampled = high_sym if len(high_sym) <= max_high_eval else random.sample(high_sym, max_high_eval)

            max_med_eval = 100
            med_sym_sampled = random.sample(med_sym, min(len(med_sym), max_med_eval))

            max_low_eval = 50
            low_sym_sampled = random.sample(low_sym, min(len(low_sym), max_low_eval))

            print(
                f"   Evaluating Physics for Sampled: High={len(high_sym_sampled)}, Med={len(med_sym_sampled)}, Low={len(low_sym_sampled)}")

            data_points = []

            def process_list(sub_list, category_name):
                for idx, res in enumerate(sub_list):
                    stab = res['stabilizer_size']
                    sym_index_rho = stab / gamma_size

                    temp_cif = os.path.join(output_dir, f"temp_{topology_name}_{idx}.cif")
                    E = evaluate_physical_energy(assembler, res['coloring'], graph_data, temp_cif, max_atoms=2500)

                    if E is not None:
                        data_points.append({'sym_index': sym_index_rho, 'energy': E, 'cat': category_name})
                        writer.writerow([topology_name, len(atom_choices), label, f"{sym_index_rho:.6f}", f"{E:.4f}",
                                         category_name])

                    if (idx + 1) % 10 == 0:
                        print(f"      Evaluated {idx + 1}/{len(sub_list)} in {category_name}...")

                    if os.path.exists(temp_cif): os.remove(temp_cif)
                    log_file = temp_cif.replace(".cif", ".log")
                    if os.path.exists(log_file): os.remove(log_file)

            process_list(high_sym_sampled, 'High')
            process_list(med_sym_sampled, 'Med')
            process_list(low_sym_sampled, 'Low')

            all_datasets.append({
                "label": label,
                "gamma_size": gamma_size,
                "points": data_points,
                "meta": {
                    "high_vals": sorted(list(high_vals), reverse=True),
                    "med_vals": sorted(list(med_vals), reverse=True)
                }
            })

    print(f"\n>>> Physical Data saved to {csv_path}")
    plot_results_individual(all_datasets, output_dir, samples_per_category)

def format_sym_index(value):
    return f"{value:.4f}".rstrip('0').rstrip('.')

def plot_results_individual(datasets, output_dir, samples_per_category):
    print(">>> Plotting Individual Symmetry-Stability Correlation Figures...")

    for data in datasets:
        fig, ax = plt.subplots(figsize=(7, 5))
        points = data['points']
        label = data['label']
        meta = data['meta']
        gamma_size = data['gamma_size']

        if not points:
            ax.text(0.5, 0.5, f"No Data\nAll structures exceeded max_atoms limit",
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=14, color='red', fontweight='bold')
            ax.set_title(f"Symmetry-Stability Correlation ({label})", fontsize=13, fontweight='bold')
        else:
            groups = {'Low': [], 'Med': [], 'High': []}
            for p in points: groups[p['cat']].append(p)

            for cat in ['Low', 'Med', 'High']:
                group_data = groups[cat]
                if not group_data: continue

                x = [p['sym_index'] for p in group_data]
                y = [p['energy'] for p in group_data]
                style = STYLE_MAP[cat]

                if cat == 'High':
                    rho_vals = [format_sym_index(v / gamma_size) for v in meta['high_vals']]
                    leg_label = f"Ordered: High Symmetry Index ({','.join(rho_vals)})"
                elif cat == 'Med':
                    rho_vals = [format_sym_index(v / gamma_size) for v in meta['med_vals']]
                    leg_label = f"Ordered: Med Symmetry Index ({','.join(rho_vals)})"
                else:
                    asym_rho = format_sym_index(1.0 / gamma_size)
                    leg_label = f"Asymmetric Noise (Symmetry Index = {asym_rho})"

                ax.scatter(x, y, c=style['color'], marker=style['marker'], label=leg_label,
                           s=70 if cat == 'High' else 45,
                           alpha=0.9 if cat == 'High' else 0.6,
                           edgecolors='black', linewidths=0.6)

            ax.set_title(f"Symmetry-Stability Correlation ({label})", fontsize=13, fontweight='bold')
            ax.set_xlabel(r"Symmetry Index ($\rho$)", fontsize=12)
            ax.set_ylabel("Thermodynamic Stability (CHGNet eV/atom)", fontsize=11)

            ax.set_xscale('log')

            all_x = [p['sym_index'] for p in points]
            if all_x:
                unique_x = sorted(list(set(all_x)))
                ax.set_xticks(unique_x)
                ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda val, pos: format_sym_index(val)))

            ax.grid(True, which="both", ls="--", alpha=0.3)
            ax.legend(loc='upper right', fontsize=9, framealpha=0.9)

            all_y = [p['energy'] for p in points]
            min_idx = np.argmin(all_y)
            min_p = points[min_idx]
            min_rho_str = format_sym_index(min_p['sym_index'])

            ax.annotate(f"Ground State\n($\\rho$ = {min_rho_str})",
                        xy=(min_p['sym_index'], min_p['energy']), xytext=(min_p['sym_index'], min_p['energy'] + 0.05),
                        arrowprops=dict(facecolor='black', shrink=0.05), fontsize=10, ha='center', fontweight='bold')

        plt.tight_layout()
        plot_path = os.path.join(output_dir, f"Symmetry_Stability_{label}_{samples_per_category}.png")
        plt.savefig(plot_path, dpi=300)
        plt.close(fig)
        print(f"   -> Plot saved to {plot_path}")

if __name__ == "__main__":
    EXPERIMENTS = [
        ("rht", ["6c_CdZn_1_Ch", "8c_Co_1_Ch"]),
        ("rht", ["6c_CdZn_1_Ch", "8c_Co_1_Ch", "0_mn1_new"]),
    ]

    run_physical_experiment(EXPERIMENTS, samples_per_category=30)