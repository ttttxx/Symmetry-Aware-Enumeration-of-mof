import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import seaborn as sns  # 引入 seaborn 用于更好的配色和样式

# ==========================================
# 0. 高级环境配置与专业样式
# ==========================================
sns.set_style("whitegrid", {
    "ytick.major.size": 0, "xtick.major.size": 0,
    "grid.color": "#f0f0f0",
    "grid.linestyle": "-"
})

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['lines.linewidth'] = 2.2
plt.rcParams['legend.fontsize'] = 9.5

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, '..', 'results')

# 2x2 布局
fig, axes = plt.subplots(2, 2, figsize=(16, 11))
plt.subplots_adjust(wspace=0.25, hspace=0.35)

current_palette = sns.color_palette("colorblind", 3)
colors = {
    'ffc': current_palette[0],
    'fjh': current_palette[1],
}

styles_2 = {'linestyle': '-', 'linewidth': 2.8, 'marker': 'o', 'markersize': 5, 'markevery': 0.1}
styles_3 = {'linestyle': '--', 'linewidth': 1.6, 'marker': 's', 'markersize': 4, 'markevery': 0.15}

topologies = ['ffc', 'fjh']
atom_counts = [3, 4, 5, 6]

# ==========================================
# 面板 A: 随机采样的对称覆盖率（不是 GA/随机的优化目标，仅说明对称集不会被顺带覆盖）
# ==========================================
ax1 = axes[0, 0]
for topo in topologies:
    for atoms in atom_counts:
        file_path = os.path.join(RESULTS_DIR, 'baseline', f'experiment_data_{topo}_10000_{atoms}.csv')
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            x_col = 'Sample_Index' if 'Sample_Index' in df.columns else 'evaluations'
            df[x_col] = df[x_col].apply(lambda x: x + 1 if x == 0 else x)

            style_config = styles_2 if atoms == 2 else styles_3
            ax1.plot(df[x_col], df['High_Sym_Coverage_Pct'],
                     label=f'{topo} ({atoms} atoms)', color=colors[topo], **style_config)

ax1.set_xscale('log')
ax1.set_xlim(1, 45000)
ax1.set_ylim(-5, 105)
ax1.set_xlabel('Number of Random Samples (Log Scale)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Symmetric Target Coverage (%)', fontsize=12, fontweight='bold')
ax1.set_title('A. Random Sampling: Symmetric Coverage', fontsize=14, fontweight='bold', pad=15)
ax1.yaxis.set_major_formatter(mtick.PercentFormatter())
ax1.legend(loc='upper left', fontsize=9, frameon=True, edgecolor='black', ncol=2)

# ==========================================
# 面板 B: 能量制导 GA 的对称覆盖率（覆盖率非其目标，仅说明对称集不会被顺带枚举）
# ==========================================
ax2 = axes[0, 1]
for topo in topologies:
    atoms_fill = 2
    file_path_fill = os.path.join(RESULTS_DIR, 'ga_baseline', f'experiment_data_{topo}_45000_{atoms_fill}.csv')
    if os.path.exists(file_path_fill):
        df_fill = pd.read_csv(file_path_fill)
        if 'Generation' not in df_fill.columns:
            df_fill['Generation'] = np.arange(1, len(df_fill) + 1)
        ax2.fill_between(df_fill['Generation'], 0, df_fill['High_Sym_Coverage_Pct'],
                         color=colors[topo], alpha=0.08, label='_nolegend_')

for topo in topologies:
    for atoms in atom_counts:
        file_path = os.path.join(RESULTS_DIR, 'ga_baseline', f'experiment_data_{topo}_45000_{atoms}.csv')
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            if 'Generation' not in df.columns:
                df['Generation'] = np.arange(1, len(df) + 1)

            style_config = styles_2 if atoms == 2 else styles_3
            style_config_B = style_config.copy()
            style_config_B['markersize'] = style_config['markersize'] + 0.5
            style_config_B['markevery'] = 15

            ax2.plot(df['Generation'], df['High_Sym_Coverage_Pct'],
                     label=f'{topo} ({atoms} atoms)', color=colors[topo], **style_config_B)

ax2.set_xlim(0, 150)
ax2.set_ylim(-5, 105)
ax2.set_xlabel('GA Generations', fontsize=12, fontweight='bold')
ax2.set_ylabel('Symmetric Target Coverage (%)', fontsize=12, fontweight='bold')
ax2.set_title('B. GA Coverage (energy-guided, not its objective)', fontsize=14, fontweight='bold', pad=15)
ax2.yaxis.set_major_formatter(mtick.PercentFormatter())

# ==========================================
# 面板 C: GA 在其自身目标上的收敛 —— 最优能量 vs 评估次数，并与枚举全局最小对比
# 面板 D: 同构冗余（核心公平论点）
# ==========================================
ax3 = axes[1, 0]  # Panel C: energy convergence
ax4 = axes[1, 1]  # Panel D: isomorphic redundancy

topo_case = 'ffc'
atom_case = 5
file_path_ga = os.path.join(RESULTS_DIR, 'ga_baseline', f'experiment_data_{topo_case}_45000_{atom_case}.csv')

PENALTY = 9000.0  # 组装失败惩罚阈值（>=此值视为无效，作图时剔除）

if os.path.exists(file_path_ga):
    df_ga = pd.read_csv(file_path_ga)

    if 'Total_Evaluations' not in df_ga.columns:
        df_ga['Total_Evaluations'] = np.arange(1, len(df_ga) + 1)
    if 'Generation' not in df_ga.columns:
        df_ga['Generation'] = np.arange(1, len(df_ga) + 1)

    color_cov = colors[topo_case]
    color_ref = '#d62728'

    # --- 面板 C: 最优能量收敛 + 枚举全局最小参考线 ---
    e_best = pd.to_numeric(df_ga['Best_Energy'], errors='coerce')
    e_best = e_best.where(e_best < PENALTY)  # 剔除惩罚值，避免纵轴被拉飞

    ax3.plot(df_ga['Total_Evaluations'], e_best,
             color=color_cov, linewidth=3.0, label='GA best energy', alpha=0.95)

    # 参考线：枚举集合的全局最小能量，由 ga_baseline.py 写入同一 CSV 的 GT_Min_Energy 常量列
    gmin = None
    gstab = None
    if 'GT_Min_Energy' in df_ga.columns:
        col = pd.to_numeric(df_ga['GT_Min_Energy'], errors='coerce').dropna()
        if len(col) > 0 and col.iloc[0] < PENALTY:
            gmin = float(col.iloc[0])
        if 'GT_Min_Stab' in df_ga.columns:
            scol = pd.to_numeric(df_ga['GT_Min_Stab'], errors='coerce').dropna()
            if len(scol) > 0:
                gstab = int(scol.iloc[0])

    if gmin is not None:
        ref_label = 'Enumerated global minimum'
        if gstab is not None:
            ref_label += f' (|Stab|={gstab})'
        ax3.axhline(gmin, color=color_ref, linestyle='--', linewidth=2.0, label=ref_label)
        # 阴影标出 GA 与枚举最优之间的能量差（未达/浪费部分）
        ax3.fill_between(df_ga['Total_Evaluations'], gmin, e_best,
                         where=(e_best > gmin), color=color_ref, alpha=0.10, label='_nolegend_')
    else:
        print("[Panel C] GT_Min_Energy column not found or empty in the CSV; reference line omitted. "
              "Run the updated ga_baseline.py (it writes GT_Min_Energy into the same CSV).")

    # 设定纵轴范围
    valid = e_best.dropna()
    if len(valid) > 0:
        lo = valid.min() if gmin is None else min(valid.min(), gmin)
        hi = valid.max()
        pad = max(0.02 * abs(hi - lo), 0.02)
        ax3.set_ylim(lo - pad, hi + pad)

    ax3.set_xlim(0, df_ga['Total_Evaluations'].max())
    ax3.set_xlabel('Total Evaluations (Computational Cost)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Best Energy (eV/atom)', fontsize=12, fontweight='bold')
    ax3.set_title(f'C. Energy Convergence vs. Enumerated Optimum ({topo_case}, {atom_case})',
                  fontsize=14, fontweight='bold', pad=15)
    ax3.legend(loc='upper right', frameon=True, edgecolor='black')

    # --- 面板 D: 同构冗余率 ---
    ax4.plot(df_ga['Generation'], df_ga['Redundancy_Rate'] * 100,
             color=color_ref, linewidth=2.5, linestyle='-.', label='Isomorphic Redundancy')
    ax4.fill_between(df_ga['Generation'], 0, df_ga['Redundancy_Rate'] * 100,
                     color=color_ref, alpha=0.18, label='_nolegend_')

    ax4.set_xlabel('Search Iterations (GA Generations)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Computational Waste / Redundancy (%)', color=color_ref, fontsize=12, fontweight='bold')
    ax4.set_xlim(0, df_ga['Generation'].max())
    ax4.set_ylim(-5, 105)
    ax4.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax4.set_title(f'D. Isomorphic Redundancy ({topo_case}, {atom_case})', fontsize=14, fontweight='bold', pad=15)

    redundant_area_max_y = df_ga['Redundancy_Rate'].max() * 100
    final_gen = df_ga['Generation'].iloc[-1]
    ax4.text(final_gen * 0.75, redundant_area_max_y * 0.82, "Waste Area",
             fontsize=11, fontweight='bold', color=color_ref,
             alpha=0.6, rotation=15, ha='center', style='italic')
    ax4.legend(loc='lower right', frameon=True, edgecolor='black')
else:
    print(f"[Panels C/D] GA data file not found: {file_path_ga}")

# ==========================================
# 保存矢量 PDF 图片
# ==========================================
plt.tight_layout()
output_path = os.path.join(RESULTS_DIR, 'enhanced_baseline_4panels.pdf')
plt.savefig(output_path, dpi=600, bbox_inches='tight')
print(f"Refined 4-panel Scientific Plot saved successfully to: {output_path}")