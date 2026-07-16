import os
import sys
import json
import csv
import subprocess
import shutil
import argparse
import matplotlib.pyplot as plt
from pymatgen.core import Structure

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(project_root)

try:
    from main import load_graph
    from assembler import StructureAssembler
except ImportError as e:
    print(f"[Error] Failed to import project modules: {e}")
    sys.exit(1)

results_root = os.path.join(project_root, "results")
HARTREE_TO_EV = 27.2114079527

def get_data_path(topology_name):
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")

def check_xtb_installed():
    if shutil.which("xtb") is None:
        print("\n[Critical Error] 'xtb' command not found in PATH.")
        return False
    return True

def cleanup_files(cif_path, xyz_path, work_dir):
    files_to_remove = [cif_path, xyz_path, xyz_path + ".out"]
    xtb_temps = ["xtbopt.xyz", "xtbopt.log", "charges", "wbo", "xtbrestart", "xtbtopo.mol", ".xtboptok"]
    for f in xtb_temps:
        f_path = os.path.join(work_dir, f)
        if os.path.exists(f_path):
            try:
                os.remove(f_path)
            except OSError:
                pass
    for f_path in files_to_remove:
        if os.path.exists(f_path):
            try:
                os.remove(f_path)
            except OSError:
                pass

def write_safe_xyz(cif_path, xyz_path):
    struct = Structure.from_file(cif_path)
    with open(xyz_path, 'w') as f:
        f.write(f"{len(struct)}\n")
        f.write("MOF_Structure_Cluster_For_xTB\n")
        for site in struct:
            f.write(f"{site.specie.symbol}  {site.x:.8f}  {site.y:.8f}  {site.z:.8f}\n")
        f.write("\n\n\n")

def run_xtb_calculation(xyz_path, work_dir):
    cmd = ["xtb", os.path.basename(xyz_path), "--gfn", "2", "--opt", "--uhf", "0"]
    output_log = xyz_path + ".out"
    try:
        with open(output_log, 'w') as f_out:
            subprocess.run(cmd, cwd=work_dir, stdout=f_out, stderr=subprocess.STDOUT, check=True)
        energy_hartree = None
        with open(output_log, 'r') as f:
            for line in f:
                if "TOTAL ENERGY" in line and "H" in line:
                    parts = line.split()
                    for p in parts:
                        try:
                            val = float(p)
                            if val < 0: energy_hartree = val
                        except ValueError:
                            continue
        return energy_hartree * HARTREE_TO_EV if energy_hartree else None
    except subprocess.CalledProcessError:
        print(f"   [Error] xTB command failed for {os.path.basename(xyz_path)}")
        return None

def calculate_structure_energy(struct_info, graph_data, assembler, work_dir, topology_name):
    coloring = struct_info.get('coloring', [])
    struct_id = struct_info.get('id', 'unknown')
    label = f"struct_{struct_id}"

    if not coloring: return None

    structure_map = {str(node_idx): str(val) for node_idx, val in enumerate(coloring)}

    temp_cif_name = f"{topology_name}_{label}_temp.cif"
    temp_xyz_name = f"XYZ_{label}.xyz"
    temp_cif_path = os.path.join(work_dir, temp_cif_name)
    temp_xyz_path = os.path.join(work_dir, temp_xyz_name)

    energy_ev = None
    success = False

    try:
        assembler.assemble_and_save(temp_cif_path, structure_map, graph_data)
        if not os.path.exists(temp_cif_path): return None

        write_safe_xyz(temp_cif_path, temp_xyz_path)
        print(f"   -> Running xTB optimization for ID {struct_id} ...")

        energy_ev = run_xtb_calculation(temp_xyz_path, work_dir)
        if energy_ev is not None:
            success = True
    except Exception as e:
        print(f"[Error] Failed processing ID {struct_id}: {e}")
    finally:
        if success: cleanup_files(temp_cif_path, temp_xyz_path, work_dir)

    return energy_ev

def run_analysis_single_file(json_path):
    if not check_xtb_installed(): return
    if not os.path.exists(json_path):
        print(f"[Error] Target JSON file not found: {json_path}")
        return

    print(f"[-] Loading JSON file: {json_path}")
    with open(json_path, 'r') as f:
        data = json.load(f)

    topology_name = data.get('meta', {}).get('topology', 'sod')
    structures_list = data.get('structures', [])

    if not structures_list:
        print("[Error] No 'structures' array found in the JSON file.")
        return

    topology_file = get_data_path(topology_name)
    xtb_work_dir = os.path.join(results_root, f"{topology_name}_xtb_workdir")
    if not os.path.exists(xtb_work_dir): os.makedirs(xtb_work_dir)

    print(f"[-] Loading Topology Graph from: {topology_file}")
    graph_data = load_graph(topology_file)

    print("[-] Initializing Assembler...")
    assembler = StructureAssembler(output_dir=xtb_work_dir)

    print(f"[-] Starting GFN-xTB Calculations for {len(structures_list)} structures...")

    max_sym = max(s.get('stabilizer_size', 1) for s in structures_list) if structures_list else 1
    data_points = []

    for i, struct_info in enumerate(structures_list):
        energy = calculate_structure_energy(struct_info, graph_data, assembler, xtb_work_dir, topology_name)
        if energy is None: continue

        sym_order = struct_info.get('stabilizer_size', 1)
        label = f"ID_{struct_info.get('id')}"

        if sym_order == max_sym:
            stype = "High"
        elif sym_order <= 2:
            stype = "Low"
        else:
            stype = "Med"

        print(f"    [{i + 1}/{len(structures_list)}] {label}: {energy:.4f} eV (Sym: {sym_order})")

        data_points.append({
            "filename": os.path.basename(json_path),
            "label": label,
            "sym_order": sym_order,
            "energy": energy,
            "type": stype
        })

    if not data_points:
        print("[Error] No valid energies calculated.")
        return

    data_points.sort(key=lambda x: x['energy'])
    base_name = os.path.splitext(os.path.basename(json_path))[0]

    csv_path = os.path.join(results_root, f"{base_name}_xtb_analysis.csv")
    save_csv(data_points, csv_path)

    plot_scatter(data_points, topology_name, base_name)
    print("[-] All Done.")

def save_csv(data, filepath):
    headers = ["Filename", "Label", "Symmetry_Order", "Energy_eV", "Type"]
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in data:
            writer.writerow([row['filename'], row['label'], row['sym_order'], row['energy'], row['type']])
    print(f"[Output] Data saved to CSV: {filepath}")

def plot_scatter(data, topology_name, base_name):
    plt.figure(figsize=(10, 7))
    style_map = {
        'High': {'color': '#d62728', 'marker': 'o', 'label': 'Highest Symmetry'},
        'Med': {'color': '#1f77b4', 'marker': '^', 'label': 'Medium Symmetry'},
        'Low': {'color': '#7f7f7f', 'marker': 'x', 'label': 'Low/No Symmetry'}
    }

    groups = {'High': [], 'Med': [], 'Low': []}
    for p in data: groups[p['type']].append(p)

    if not data: return
    min_e = min(p['energy'] for p in data)

    for g_name in ['High', 'Med', 'Low']:
        points = groups[g_name]
        if not points: continue
        x = [p['sym_order'] for p in points]
        y = [p['energy'] - min_e for p in points]
        style = style_map[g_name]

        kwargs = {'c': style['color'], 'marker': style['marker'], 'label': style['label'], 's': 120, 'alpha': 0.8}
        if style['marker'] != 'x':
            kwargs['edgecolors'] = 'k'
        else:
            kwargs['linewidths'] = 2

        plt.scatter(x, y, **kwargs)

    plt.title(f'Symmetry vs Stability: {topology_name.upper()} ({base_name})', fontsize=14, fontweight='bold')
    plt.xlabel('Symmetry Order (Stabilizer Size)', fontsize=14)
    plt.ylabel('Relative Energy (eV)', fontsize=14)
    plt.xscale('log')
    plt.grid(True, linestyle='--', alpha=0.4, which="both")
    plt.legend(fontsize=12)
    plt.tight_layout()

    output_img = os.path.join(results_root, f"{base_name}_xtb_analysis.png")
    plt.savefig(output_img, dpi=300)
    print(f"[Output] Plot saved to PNG: {output_img}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run xTB analysis on a specific JSON file containing multiple structures.")
    parser.add_argument("-f", "--file", type=str, required=True,
                        help="Path to the JSON file to analyze (e.g., results/sod_0_mn1_new6c_Al_1_Ch_20260108_193244.json)")

    args = parser.parse_args()

    os.chdir(project_root)

    run_analysis_single_file(args.file)
