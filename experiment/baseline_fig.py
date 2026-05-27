import os
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

# 修改核心1：将 1x3 布局改为 2x2 布局以容纳拆分后的图
fig, axes = plt.subplots(2, 2, figsize=(16, 11))
plt.subplots_adjust(wspace=0.25, hspace=0.35)

current_palette = sns.color_palette("colorblind", 3)
colors = {
    'sod': current_palette[0],
    'lil': current_palette[1],
    'rht': current_palette[2]
}

styles_2 = {'linestyle': '-', 'linewidth': 2.8, 'marker': 'o', 'markersize': 5, 'markevery': 0.1}
styles_3 = {'linestyle': '--', 'linewidth': 1.6, 'marker': 's', 'markersize': 4, 'markevery': 0.15}

topologies = ['sod', 'lil', 'rht']
atom_counts = [2, 3]

# ==========================================
# 面板 A: 随机搜索的“空间稀释”
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
ax1.set_title('A. Stochastic Search Complexity', fontsize=14, fontweight='bold', pad=15)
ax1.yaxis.set_major_formatter(mtick.PercentFormatter())
ax1.legend(loc='upper left', fontsize=9, frameon=True, edgecolor='black', ncol=2)

# ==========================================
# 面板 B: 遗传算法的覆盖率瓶颈
# ==========================================
ax2 = axes[0, 1]
for topo in topologies:
    atoms_fill = 2
    file_path_fill = os.path.join(RESULTS_DIR, 'ga_baseline', f'experiment_data_{topo}_45000_{atoms_fill}.csv')
    if os.path.exists(file_path_fill):
        df_fill = pd.read_csv(file_path_fill)
        # 兼容处理 Generation
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
ax2.set_title('B. GA Performance Ceiling', fontsize=14, fontweight='bold', pad=15)
ax2.yaxis.set_major_formatter(mtick.PercentFormatter())

# ==========================================
# 面板 C & 面板 D: 拆分“同构陷阱”的双轴图
# ==========================================
ax3 = axes[1, 0]  # 新的面板 C
ax4 = axes[1, 1]  # 新的面板 D

topo_case = 'rht'  # 保持展示瓶颈最严重的 rht
atom_case = 2
file_path_ga = os.path.join(RESULTS_DIR, 'ga_baseline', f'experiment_data_{topo_case}_45000_{atom_case}.csv')

if os.path.exists(file_path_ga):
    df_ga = pd.read_csv(file_path_ga)

    # 动态数据兼容：如果你的旧CSV没这两个字段，我们向下兼容计算
    if 'Total_Evaluations' not in df_ga.columns:
        pop_size = len(df_ga)  # fallback
        df_ga['Total_Evaluations'] = np.arange(1, len(df_ga) + 1) * pop_size
    if 'Generation' not in df_ga.columns:
        df_ga['Generation'] = np.arange(1, len(df_ga) + 1)

    color_cov = colors[topo_case]
    color_red = '#d62728'

    # --- 绘制图 C: 计算代价与性能天花板 ---
    ax3.plot(df_ga['Total_Evaluations'], df_ga['High_Sym_Coverage_Pct'],
             color=color_cov, linewidth=3.2, label='Target Coverage', alpha=0.9)
    ax3.set_xlabel('Total Evaluations (Computational Cost)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Symmetric Target Coverage (%)', color=color_cov, fontsize=12, fontweight='bold')
    ax3.set_xlim(0, max(df_ga['Total_Evaluations']))
    ax3.set_ylim(-5, 105)
    ax3.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax3.set_title(f'C. Evaluation Efficiency Ceiling ({topo_case}, {atom_case})', fontsize=14, fontweight='bold',
                  pad=15)
    ax3.legend(loc='lower right', frameon=True, edgecolor='black')

    # --- 绘制图 D: 搜索迭代次数与冗余率 (真正的搜索迭代分析) ---
    ax4.plot(df_ga['Generation'], df_ga['Redundancy_Rate'] * 100,
             color=color_red, linewidth=2.5, linestyle='-.', label='Isomorphic Redundancy')
    # 红色面积填充表示浪费计算量
    ax4.fill_between(df_ga['Generation'], 0, df_ga['Redundancy_Rate'] * 100,
                     color=color_red, alpha=0.18, label='_nolegend_')

    ax4.set_xlabel('Search Iterations (GA Generations)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Computational Waste / Redundancy (%)', color=color_red, fontsize=12, fontweight='bold')
    ax4.set_xlim(0, max(df_ga['Generation']))
    ax4.set_ylim(-5, 105)
    ax4.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax4.set_title(f'D. The "Isomorphism Trap" ({topo_case}, {atom_case})', fontsize=14, fontweight='bold', pad=15)

    # 在图 D 中添加 Waste Area 标注
    redundant_area_max_y = df_ga['Redundancy_Rate'].max() * 100
    final_gen = df_ga['Generation'].iloc[-1]
    ax4.text(final_gen * 0.75, redundant_area_max_y * 0.82, "Waste Area",
             fontsize=11, fontweight='bold', color=color_red,
             alpha=0.6, rotation=15, ha='center', style='italic')
    ax4.legend(loc='lower right', frameon=True, edgecolor='black')

# ==========================================
# 保存矢量 PDF 图片
# ==========================================
plt.tight_layout()
output_path = os.path.join(RESULTS_DIR, 'enhanced_baseline_4panels.pdf')
plt.savefig(output_path, dpi=600, bbox_inches='tight')
print(f"Refined 4-panel Scientific Plot saved successfully to: {output_path}")