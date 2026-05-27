### Input
``data/mof/tbo.cif``

### Step 1 (Preprocessing)

run``python utils/parse_cif.py data/mof/tbo.cif ``

内部调用 external_lib 的手动解析逻辑读取 CIF。

内部调用 external_lib.utils.periodic_structure_analysis 计算原子间连接（Build Graph）。

Output: 生成 data/tbo/tbo_graph.json（包含节点、边、分数坐标、师姐逻辑下的 Symmetry 标记）。

### Step 2 (Analysis)
`` run mof/main.py ``

main.py 读取 tbo_graph.json。

SymmetryEngine 使用师姐的 PBC 算法验证距离一致性。

### Step 3 实验
结构稳定性
先 ``run structure_selection.py``得到不同稳定性的几组json
再``run run_analysis_plot.py``进行分析

baseline
直接``run baseline.py``得到baseline的冗余率balabala，直接用于对照
