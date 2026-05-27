import os
import sys
import random
import time
import csv
import numpy as np

# [配置] 解决 macOS 下 matplotlib 报错并设置为后台绘图
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 路径修正：确保能导入上级目录的模块
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(current_dir)

try:
    from main import load_graph
    from symmetry import SymmetryEngine
    from enumerator import MOFEnumerator
except ImportError:
    pass


def get_data_path(topology_name):
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")


def run_sampling_experiment_improved(topology_name, atom_choices, num_samples=5000):
    print(f"\n>>> Running Baseline Experiment on [{topology_name}]")
    print(f"    Atom Choices ({len(atom_choices)} types): {atom_choices}")

    graph_path = get_data_path(topology_name)
    if not os.path.exists(graph_path):
        print(f"Error: Graph file not found at {graph_path}")
        return [], [], []

    graph_data = load_graph(graph_path)

    # --- 1. 获取 Ground Truth (理论全集) ---
    print("1. Computing Ground Truth (Burnside)...")
    sym = SymmetryEngine(graph_data, verbose=False)
    perms = sym.compute_automorphism_group()

    enum = MOFEnumerator(perms, graph_data['n_nodes'], verbose=False)
    # 获取所有唯一结构的集合（作为参考系）
    all_unique_structures = enum.enumerate(atom_choices)

    # 建立查找表：Canonical -> 是否高对称
    # [定义] 高对称性：stabilizer_size > 1 (非对称结构 stab=1)
    unique_orbits_map = {}
    high_sym_canonicals = set()

    for res in all_unique_structures:
        canon = tuple(res['coloring'])
        unique_orbits_map[canon] = res
        if res['stabilizer_size'] > 1:
            high_sym_canonicals.add(canon)

    total_unique_count = len(all_unique_structures)
    total_high_sym_targets = len(high_sym_canonicals)  # 总高对称目标数

    print(f"   -> Ground Truth: {total_unique_count} unique orbits.")
    print(f"   -> High-Sym Targets (Stab > 1): {total_high_sym_targets}")

    # --- 2. 随机采样与冗余分析 ---
    print(f"2. Random Sampling {num_samples} structures...")

    found_canonicals = set()  # 记录已发现的唯一结构
    found_history = []  # 历史数据：发现唯一数量
    redundancy_history = []  # 历史数据：冗余率

    # [修改] 记录高对称覆盖率 (%) 而非数量
    high_sym_coverage_history = []

    nodes = graph_data['n_nodes']
    atoms_pool = sorted(atom_choices)

    start_time = time.time()

    for k in range(1, num_samples + 1):
        # A. 随机生成
        random_coloring = tuple(random.choices(atoms_pool, k=nodes))

        # B. 正则化 (Canonicalization)
        min_coloring = random_coloring
        for p in perms:
            mapped = tuple(random_coloring[i] for i in p)
            if mapped < min_coloring:
                min_coloring = mapped

        # C. 记录
        found_canonicals.add(min_coloring)

        # D. 统计
        unique_so_far = len(found_canonicals)
        found_history.append(unique_so_far)

        # 实时冗余率
        current_redundancy = 1.0 - (unique_so_far / k)
        redundancy_history.append(current_redundancy)

        # [修改] 实时高对称覆盖率 (Rate)
        current_high_sym_count = sum(1 for c in found_canonicals if c in high_sym_canonicals)

        if total_high_sym_targets > 0:
            current_coverage_pct = (current_high_sym_count / total_high_sym_targets) * 100.0
        else:
            current_coverage_pct = 0.0

        high_sym_coverage_history.append(current_coverage_pct)

    elapsed = time.time() - start_time

    # --- 3. 结果汇总 ---
    final_high_sym_pct = high_sym_coverage_history[-1]

    global_cov = len(found_canonicals) / total_unique_count if total_unique_count > 0 else 0
    final_red = redundancy_history[-1]

    print(f"   -> Time: {elapsed:.2f}s")
    print(f"   -> Metric 1 (Global Coverage):   {len(found_canonicals)}/{total_unique_count} ({global_cov * 100:.2f}%)")
    print(f"   -> Metric 2 (High-Sym Coverage): {final_high_sym_pct:.2f}%")
    print(f"   -> Metric 3 (Redundancy):        {final_red * 100:.2f}%")

    # 返回所有历史数据以便保存
    return found_history, redundancy_history, high_sym_coverage_history


if __name__ == "__main__":
    # ================= 配置区 =================
    TOPOLOGY = "sod"  # 拓扑名称
    SAMPLES = 10000  # 采样次数
    # 原子组合
    ATOMS = ["6c_CdZn_1_Ch", "8c_Co_1_Ch", "0_mn1_new"]
    # ATOMS = ["6c_CdZn_1_Ch", "8c_Co_1_Ch"]
    # =========================================


    atom_count = len(ATOMS)

    # 运行实验
    history_unique, history_red, history_high_cov = run_sampling_experiment_improved(TOPOLOGY, ATOMS,
                                                                                     num_samples=SAMPLES)

    if history_unique:
        # 1. 准备目录
        output_dir = os.path.join(project_root, "results", "baseline")
        os.makedirs(output_dir, exist_ok=True)

        # 2. 生成文件名: experiment_data_{TOPOLOGY}_{SAMPLES}_{ATOMS_COUNT}.csv
        filename_base = f"experiment_data_{TOPOLOGY}_{SAMPLES}_{atom_count}"
        csv_path = os.path.join(output_dir, f"{filename_base}.csv")

        print(f"\n>>> Saving data to {csv_path}...")

        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            # [修改] CSV Header 更新为覆盖率百分比
            writer.writerow(["Sample_Index", "Unique_Found", "Redundancy_Rate", "High_Sym_Coverage_Pct"])
            for i in range(SAMPLES):
                writer.writerow([
                    i + 1,
                    history_unique[i],
                    history_red[i],
                    history_high_cov[i]
                ])

        # 3. 绘图并保存
        print(">>> Plotting...")
        plt.figure(figsize=(4, 4))

        # 子图3: 高对称覆盖率 (High-Symmetry Target Coverage)
        plt.subplot(1, 1, 1)
        plt.plot(history_high_cov, color='g', label='Stochastic Search')
        plt.xlabel('Number of Random Samples')
        # [修改] Y轴标签改为百分比
        plt.ylabel('Symmetric Target Coverage (%)')
        plt.title(f'{TOPOLOGY}, {atom_count} atoms')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        img_path = os.path.join(output_dir, f"{filename_base}.png")
        plt.savefig(img_path, dpi=300)
        print(f"   -> Plot Saved to {img_path}")