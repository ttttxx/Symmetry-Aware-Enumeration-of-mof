import os
import json
import glob
import pandas as pd
import matplotlib.pyplot as plt


def main():
    # ================= 配置区 =================
    DATA_PREFIX = "sod"  # 请与你 main.py 中的设置保持一致
    J_PENALTY = 1.0  # Ising 模型的同种原子排斥惩罚系数 (J > 0)
    # =========================================

    print(f"--- Generating Energy Landscape for {DATA_PREFIX} ---")

    # 1. 读取拓扑图以获取边 (Edges) 信息，用于计算相互作用
    graph_path = os.path.join("data", DATA_PREFIX, f"{DATA_PREFIX}_graph.json")
    if not os.path.exists(graph_path):
        raise FileNotFoundError(f"找不到拓扑图文件: {graph_path}")

    with open(graph_path, 'r') as f:
        graph_data = json.load(f)
    edges = graph_data["edges"]

    # 2. 自动寻找 results 目录下最新的 json 结果文件
    result_files = glob.glob(os.path.join("results", f"{DATA_PREFIX}_*.json"))
    if not result_files:
        raise FileNotFoundError("没有找到结果文件，请确保已经成功运行了 main.py！")

    latest_result_file = max(result_files, key=os.path.getctime)
    print(f"Loading enumeration results from: {latest_result_file}")

    with open(latest_result_file, 'r') as f:
        result_data = json.load(f)

    # 提取全群大小 |Gamma| 和构型列表
    group_size = result_data["meta"]["group_size"]
    structures = result_data["structures"]

    # 3. 遍历构型，计算对称性指数 rho(c) 和 Ising 能量 E(c)
    data_records = []
    for struct in structures:
        coloring = struct["coloring"]
        stab_size = struct["stabilizer_size"]

        # 计算对称性指数 rho(c) = |Stab(c)| / |Gamma|
        rho = stab_size / group_size

        # 计算 Ising-type Lattice Energy
        # 遍历所有边，如果两端原子相同，则加上排斥惩罚 J
        energy = 0.0
        for edge in edges:
            u = edge["source"]
            v = edge["target"]
            if coloring[u] == coloring[v]:
                energy += J_PENALTY

        data_records.append({
            "id": struct["id"],
            "rho": rho,
            "energy": energy,
            "stabilizer_size": stab_size
        })

    # 4. 保存为 CSV，供后续查阅或过滤
    df = pd.DataFrame(data_records)
    csv_path = os.path.join("results", f"{DATA_PREFIX}_energy_scatter.csv")
    df.to_csv(csv_path, index=False)
    print(f"Data saved to {csv_path}")

    # ================= 5. 绘制高颜值的 SCI 矢量散点图 =================
    plt.figure(figsize=(10, 5))

    # 区分不对称极限 (Asymmetric Limit)
    asymmetric_limit = 1.0 / group_size

    # 将数据分为两组：不对称噪声 (灰色) 和 对称极小值 (彩色)
    mask_asym = df['rho'] <= asymmetric_limit + 1e-6
    mask_sym = df['rho'] > asymmetric_limit + 1e-6

    # 绘制不对称点 (底层，半透明灰色)
    plt.scatter(df[mask_asym]['rho'], df[mask_asym]['energy'],
                c='#B0BEC5', alpha=0.4, s=20, label='Asymmetric Limit (Filtered)')

    # 绘制对称点 (顶层，按能量高低赋予冷暖色调)
    scatter = plt.scatter(df[mask_sym]['rho'], df[mask_sym]['energy'],
                          c=df[mask_sym]['energy'], cmap='viridis_r',
                          edgecolor='black', linewidth=0.5, s=70, alpha=0.9,
                          label='Symmetric Candidates', zorder=5)

    # 设置图表样式
    plt.xlabel(r'Symmetry Index $\rho(c)$', fontsize=14, fontweight='bold')
    plt.ylabel('Lattice Energy $E(c)$', fontsize=14, fontweight='bold')
    plt.title('Symmetry-Stability Correlation', fontsize=15, fontweight='bold')

    # 坐标轴与网格优化
    plt.grid(True, linestyle='--', alpha=0.5, zorder=0)
    plt.colorbar(scatter, label='Repulsion Energy (Lower is more stable)')
    plt.legend(loc='upper right', framealpha=0.9)

    # 保存为透明底的 SVG 矢量图
    svg_path = os.path.join("results", f"{DATA_PREFIX}_energy_scatter.svg")
    plt.savefig(svg_path, transparent=True, bbox_inches='tight', format='svg', dpi=300)
    print(f"Plot saved to {svg_path}")

    # 如果你想直接看一眼图，取消下面这行的注释
    # plt.show()


if __name__ == "__main__":
    main()