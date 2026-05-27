import os
import sys
import json
import csv
import matplotlib.pyplot as plt
import numpy as np

# --- 1. 路径与环境设置 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(project_root)

# 导入 main 中的 load_graph 以获取拓扑的边(Edges)信息
try:
    from main import load_graph
except ImportError:
    # 如果导入失败，提供一个简单的本地加载函数
    def load_graph(path):
        with open(path, 'r') as f:
            return json.load(f)

results_root = os.path.join(project_root, "results", "structure_selection")

# --- 2. 物理模型参数 (可根据化学直觉调整) ---
# 假设系统倾向于形成异种原子键 (A-B)，因为这通常能最大化库仑相互作用或应力释放
# 能量单位: eV (相对值)
BOND_ENERGIES = {
    "same_type": 0.10,  # A-A 或 B-B 键 (同种排斥/弱键，能量较高)
    "diff_type": -0.20  # A-B 键 (异种吸引/强键，能量较低)
}


def get_data_path(topology_name):
    # 假设数据在 mof/data/sod/sod_graph.json
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")


def calculate_ising_energy(structure_data, graph_data):
    """
    【基于伊辛模型的真实能量计算】
    计算晶体中所有相邻节点的键能总和。
    兼容 List [u, v] 和 Dict {'source': u, 'target': v} 两种格式。
    """
    coloring = structure_data['coloring']  # 例如 [0, 1, 0, 1...]

    # 获取边列表或邻接表
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

    # Case A: 处理显式的边列表 (Edges List)
    if edges_raw is not None:
        for edge in edges_raw:
            u, v = None, None

            # --- [智能兼容修复] ---
            if isinstance(edge, dict):
                # 尝试读取常见的字典键名
                # 常见格式: {'source': 0, 'target': 1} 或 {'u': 0, 'v': 1}
                u = edge.get('source', edge.get('u'))
                v = edge.get('target', edge.get('v'))
                # 如果还是None，尝试 'node1', 'node2'
                if u is None: u = edge.get('node1')
                if v is None: v = edge.get('node2')
            elif isinstance(edge, (list, tuple)):
                # 列表格式 [0, 1]
                if len(edge) >= 2:
                    u = edge[0]
                    v = edge[1]

            # 如果没读到有效的 u, v，跳过
            if u is None or v is None:
                continue

            # 安全检查：防止索引越界
            if u >= len(coloring) or v >= len(coloring):
                continue

            type_u = coloring[u]
            type_v = coloring[v]

            if type_u == type_v:
                total_energy += BOND_ENERGIES['same_type']
            else:
                total_energy += BOND_ENERGIES['diff_type']
            valid_edges_count += 1

        # 如果读取到了边，直接返回；如果没读到（可能是格式太怪），尝试 fallback 到邻接表
        if valid_edges_count > 0:
            return round(total_energy, 4)

    # Case B: 处理邻接表 (Adjacency List)
    # 格式通常是 [[neighbor_idxs...], [neighbor_idxs...]]
    if adj_list is not None:
        for u, neighbors in enumerate(adj_list):
            type_u = coloring[u]
            for val in neighbors:
                # 邻接表里的项也可能是 (v, dist) 或纯 v 或 dict
                v = None
                if isinstance(val, (int, np.integer)):
                    v = int(val)
                elif isinstance(val, (list, tuple)):
                    v = int(val[0])
                elif isinstance(val, dict):
                    v = val.get('id', val.get('node', None))

                if v is None: continue

                # 避免重复计算无向边 (只计算 u < v)
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
    # 1. 准备路径
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

    # 2. 遍历结构并计算
    print(f"[-] Calculating energy using Bond Valence Model (Ising)...")
    for filename in json_files:
        filepath = os.path.join(target_dir, filename)

        with open(filepath, 'r') as f:
            struct_data = json.load(f)

        # === 核心调用 ===
        energy = calculate_ising_energy(struct_data, graph_data)

        # 将计算结果更新到内存对象
        struct_data['properties']['energy_ev'] = energy

        # 判断类型 (High/Med/Low)
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

    # 3. 按能量排序方便查看
    data_points.sort(key=lambda x: x['energy'])

    # 4. 保存为 CSV
    csv_path = os.path.join(results_root, f"{topology_name}_analysis_data.csv")
    save_csv(data_points, csv_path)

    # 5. 绘图
    plot_scatter(data_points, topology_name)
    print("Done.")


def save_csv(data, filepath):
    """将分析结果保存为标准 CSV 文件"""
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

    # 定义样式
    style_map = {
        'High': {'color': '#d62728', 'marker': 'o', 'label': 'High Symmetry (Top-10)'},  # Red
        'Med': {'color': '#1f77b4', 'marker': '^', 'label': 'Medium Symmetry'},  # Blue
        'Low': {'color': '#7f7f7f', 'marker': 'x', 'label': 'Asymmetric (Low-10)'}  # Grey
    }

    # 分组绘制
    groups = {'High': [], 'Med': [], 'Low': []}
    for p in data:
        groups[p['type']].append(p)

    for g_name in ['High', 'Med', 'Low']:
        points = groups[g_name]
        if not points: continue

        x = [p['sym_order'] for p in points]
        y = [p['energy'] for p in points]
        style = style_map[g_name]

        # 修复 UserWarning: 对于 marker='x' (unfilled)，不要设置 edgecolors
        if style['marker'] == 'x':
            plt.scatter(x, y, c=style['color'], marker=style['marker'],
                        label=style['label'], s=120, alpha=0.8, linewidths=2)
        else:
            plt.scatter(x, y, c=style['color'], marker=style['marker'],
                        label=style['label'], s=120, alpha=0.8, edgecolors='k')

    plt.title(f'Symmetry vs Stability: {topology_name.upper()} Topology', fontsize=16, fontweight='bold')
    plt.xlabel('Symmetry Order (Size of Stabilizer Group)', fontsize=14)
    plt.ylabel('Configurational Energy (eV) [Ising Model]', fontsize=14)

    # 使用对数坐标轴
    plt.xscale('log')

    plt.grid(True, linestyle='--', alpha=0.4, which="both")
    plt.legend(fontsize=12)
    plt.tight_layout()

    output_img = os.path.join(results_root, f"{topology_name}_analysis.png")
    plt.savefig(output_img, dpi=300)
    print(f"[Output] Plot saved to PNG: {output_img}")


if __name__ == "__main__":
    run_analysis("sod")