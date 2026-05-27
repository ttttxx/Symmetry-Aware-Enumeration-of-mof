import os
import sys
import json
import numpy as np

# --- 1. 路径修正与环境设置 ---
# 获取当前脚本所在目录 (mof/experiment)
current_dir = os.path.dirname(os.path.abspath(__file__))
# 获取项目根目录 (mof)
project_root = os.path.dirname(current_dir)
# 将根目录加入路径，以便导入 main, symmetry, enumerator
sys.path.append(project_root)

from main import load_graph
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator


def get_data_path(topology_name):
    """
    根据新的目录结构获取数据路径:
    mof/experiment/script.py -> mof/data/topology/file.json
    """
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")


def ensure_output_dir(subfolder):
    """确保输出目录存在: mof/results/experiment/"""
    # 假设 results 在 mof/results
    out_path = os.path.join(project_root, "results", "structure_selection", subfolder)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    return out_path


def select_and_save_structures(topology_name, atom_choices):
    print(f"\n>>> Running Structure Selection for Analysis: {topology_name} {atom_choices}")

    # 1. 加载与计算
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

    # 2. 排序与筛选逻辑
    # 我们使用 'stabilizer_size' 作为对称性的度量指标
    # stabilizer_size 越大 -> 对称性越高
    # stabilizer_size = 1 -> 不对称 (Asymmetric)

    # 按对称性从高到低排序
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

    # 3. 保存结果
    output_dir = ensure_output_dir(f"{topology_name}_selection")

    print(f"   -> Saving {len(selected_samples)} selected structures to {output_dir}...")

    summary_list = []

    for res in selected_samples:
        label = res['_label']
        # 构造文件名: sod_HighSym_1_CdZn_Co.json
        filename = f"{topology_name}_{label}.json"
        save_path = os.path.join(output_dir, filename)

        # 准备保存的数据
        save_data = {
            "topology": topology_name,
            "label": label,
            "symmetry_order": res['stabilizer_size'],  # 对称阶数
            "coloring": res['coloring'],  # 具体的原子排列
            "atom_choices": atom_choices,
            # 预留给后续实验填写的字段
            "properties": {
                "energy_ev": None,  # 待计算
                "stability_score": None  # 待计算
            }
        }

        with open(save_path, 'w') as f:
            json.dump(save_data, f, indent=4)

        summary_list.append({
            "Label": label,
            "Sym_Order": res['stabilizer_size'],
            "Coloring_Snippet": str(res['coloring'][:6]) + "..."
        })

    # 打印一个简报表
    print(f"\n{'-' * 60}")
    print(f"{'Label':<15} | {'Sym Order':<10} | {'Structure Snippet'}")
    print(f"{'-' * 60}")
    for item in summary_list:
        print(f"{item['Label']:<15} | {item['Sym_Order']:<10} | {item['Coloring_Snippet']}")
    print(f"{'-' * 60}")
    print(f"Done. You can now use files in {output_dir} for stability analysis.\n")


if __name__ == "__main__":
    # 配置实验
    # 建议使用稍微复杂一点的组合，以便能区分出中等对称性
    # 如果是对称性极高的拓扑（如sod）配合简单原子，可能只有 High 和 Low，没有 Middle

    # 实验 1: SOD 多元组
    select_and_save_structures("sod", ["6c_CdZn_1_Ch", "8c_Co_1_Ch"])

    # 实验 2: RHT 二元组
    # select_and_save_structures("rht", ["A", "B"])