import os
import sys
import time
import csv
import math
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# --- 路径设置 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(current_dir)
sys.path.append(project_root)

from main import load_graph
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator


def get_data_path(topology_name):
    return os.path.join(project_root, "data", topology_name, f"{topology_name}_graph.json")


def benchmark_case(topology_name, atom_choices):
    """运行单个测试用例并返回详细指标"""
    path = get_data_path(topology_name)
    if not os.path.exists(path):
        print(f"Warning: {path} not found.")
        return None

    graph_data = load_graph(path)
    n_nodes = graph_data['n_nodes']
    n_atoms = len(atom_choices)

    print(f"   Benchmarking {topology_name} ({n_atoms} types)...")

    # 1. 计时：对称性分析
    t0 = time.time()
    sym = SymmetryEngine(graph_data, verbose=False)
    perms = sym.compute_automorphism_group()
    t_sym = time.time() - t0

    # 获取轨道数 (用于计算理论缩减)
    # 简单的轨道计算：用 Identity (g=e) 跑一次并查集 或者直接从 enumerator 内部机制推断
    # 这里我们快速估算：Orb ≈ Nodes / (Size_G / Size_Stab_avg) -> 比较复杂
    # 最准确的方法是实例化 Enumerator 偷看它的 orbits
    temp_enum = MOFEnumerator(perms, n_nodes, verbose=False)
    # 在 Identity 任务下的轨道数就是全对称下的独立节点数？
    # 不，我们用最保守的估计：全群作用下的轨道数
    # MOFEnumerator._compute_orbits_by_gen(identity) 会返回全离散的轨道
    # 我们需要的是 "High Symmetry Task" 的轨道数。
    # 为了简化对比，我们计算 "Naive Space" vs "Actual Output"

    # 2. 计时：枚举
    t1 = time.time()
    results = temp_enum.enumerate(atom_choices)
    t_enum = time.time() - t1

    n_unique = len(results)

    # 3. 计算理论复杂度指标
    # Naive Space: |Atom|^|Nodes|
    naive_space_log10 = n_nodes * math.log10(n_atoms)

    return {
        "Topology": topology_name,
        "Atom_Count": n_atoms,
        "Nodes": n_nodes,
        "Naive_Log10": naive_space_log10,  # 10的多少次方
        "Unique_Found": n_unique,
        "Time_Total": t_sym + t_enum
    }


def run_performance_experiment():
    print(">>> Starting Complexity & Performance Benchmark...")

    # 定义测试集：(拓扑, 原子数量)
    # 我们通过增加原子种类来模拟复杂度提升
    test_cases = [
        ("sod", 2),
        ("sod", 3),
        ("rht", 2),
        ("rht", 3),
        ("lil", 2),
        ("lil", 3),
        # ("fcu", 2), # 如果有 fcu 数据建议加上，效果最震撼
        # ("tbo", 2)  # 如果有 tbo 数据建议加上
    ]

    # 构造原子列表
    def get_atoms(n):
        base = ["6c_CdZn_1_Ch", "8c_Co_1_Ch", "0_mn1_new", "12c_Eu_1_Ch", "4c_Ni_1_Ch"]
        return base[:n]

    results = []

    output_dir = os.path.join(project_root, "results", "performance")
    os.makedirs(output_dir, exist_ok=True)

    for topo, n_atm in test_cases:
        atoms = get_atoms(n_atm)
        data = benchmark_case(topo, atoms)
        if data:
            results.append(data)

    # 保存 CSV
    csv_path = os.path.join(output_dir, "performance_data.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f">>> Data saved to {csv_path}")
    plot_performance(results, output_dir)


def plot_performance(data, output_dir):
    print(">>> Plotting Performance Chart...")

    labels = [f"{d['Topology']}-{d['Atom_Count']}" for d in data]
    times = [d['Time_Total'] for d in data]
    naive_logs = [d['Naive_Log10'] for d in data]

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # 左轴：实际运行时间 (Bar)
    color = 'tab:blue'
    ax1.set_xlabel('Test Case (Topology - Atom Types)', fontweight='bold')
    ax1.set_ylabel('Actual Runtime (seconds)', color=color, fontweight='bold')
    bars = ax1.bar(labels, times, color=color, alpha=0.6, width=0.5, label='Actual Time')
    ax1.tick_params(axis='y', labelcolor=color)

    # 标注时间
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2., height,
                 f'{height:.2f}s',
                 ha='center', va='bottom', color='black', fontsize=9)

    # 右轴：理论搜索空间 (Line)
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Naive Search Space ($10^N$)', color=color, fontweight='bold')
    ax2.plot(labels, naive_logs, color=color, marker='o', linewidth=2, linestyle='--', label='Naive Space')
    ax2.tick_params(axis='y', labelcolor=color)

    # 标注指数
    for i, val in enumerate(naive_logs):
        ax2.text(i, val, f'$10^{{{val:.1f}}}$', ha='center', va='bottom', color='red', fontsize=9)

    plt.title('Performance vs. Combinatorial Explosion', fontsize=14, fontweight='bold')
    fig.tight_layout()

    plot_path = os.path.join(output_dir, "complexity_performance.png")
    plt.savefig(plot_path, dpi=300)
    print(f"   -> Plot saved to {plot_path}")


if __name__ == "__main__":
    run_performance_experiment()