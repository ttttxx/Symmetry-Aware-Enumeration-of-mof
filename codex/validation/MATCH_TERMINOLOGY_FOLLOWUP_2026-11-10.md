# 术语校准追查与第一轮报告更正

日期：2026-11-10。承接作者针对第一轮报告提出的四项问题。本轮只核查和提出具体改法，未修改主稿、正式数据或相邻 CO2_N2 项目的代码。

后续执行状态：作者已决定 topology recovery 不改。PLD/LCD 解析器、当前论文描述符及九组正文范围已修正，相关 4 项测试通过；模板表已移入附录 A.1，补齐十种节点和固定连接体。下文“未修改”描述的是本追查报告形成时的阶段；执行结果以 PROJECT_STATE 和 WORKLOG 为准。

## 1. 已经用对的术语不应列成待纠错项

- `node assignment` 与 `coloring`：2.1 已明确定义为同一映射的化学任务描述与数学表示，使用合适。无需强制全文只保留其中一个词，也无需再次要求作者敲定两者的定义。
- `ordered arrangement` / `ordering`：用来讨论类型的空间排布是合适的；只在指代算法输出时以 canonical representative 更直接。变化本身不构成概念混用。
- `degree-compatible`：英文准确，按图顶点度数匹配节点连接位点数；Results 相应中文“配位数相容”需改为“度数相容”。
- `stabilizer fraction`：定义和匹配分组句均准确。DFT 结果的 `symmetry index` 泛指同一个量，建议具体写 stabilizer fraction；这属于精确程度改进，并非公式或概念定义错误。
- `geometric construction`：英语主叙述已经基本统一；内部 `reconstruction` 标签和代码文件名无需为了表面统一而改名。
- `building block` 和 `building unit` 可在一般 MOF 叙述中表达相同总称，不能仅因两个词出现就定为错误；`secondary building unit` 是否覆盖所有节点才是需具体判断的化学层级问题。

第一轮最终回复把上述“现有区分应保留”与“实际需修改”放在同一待敲定表中，容易误导。后续逐句记录应明确标记：合适且保留、可选精化、概念/证据需修改。

## 2. topology recovery：实际调用链与建议文本

- `experiment/geometric_reconstruction_summary.py` 的审计使用 `StructureValidator(strict_mode=False, check_symmetry=False)` 和 `run_symmetry=False`。
- `assembler.py:792` 调用验证器，传入 `topology_graph=TG`、`placed_all=fc_placed_all`、`bonds_all=fixed_bonds`、`lattice_matrix=sc_unit_cell`。
- `StructureValidator.validate` 调用 `_check_topology`（`structure_validator.py:314`）。该函数用构建单元归属与传入键表恢复 node–linker–node 连接，将组装图和目标图转为简单图，再比较节点数、边数、度序列和图同构。
- `_check_topology` 的函数注释称“从 3D 坐标重新推导”，但函数体未直接用 `structure` 或 `lattice_matrix` 重新识别键，也未保留周期平移标记。不能依据该注释扩大检查范围。
- `fc_placed_all` 和 `fixed_bonds` 来自组装流程；此检查并非从清理后导出的 CIF 独立恢复连接。最终结构的几何短接触检查另有实现，须分别描述。

建议摘要将 `models recovered the target topology` 改为：

> the models passed a building-unit connectivity check

建议 Results 将 `recovered the target topology` 改为：

> its building-unit connectivity graph was isomorphic to the simplified target graph

可形成完整句：

> Each model yielded a sanitized, reloadable CIF and passed the short-contact check. Its building-unit connectivity graph was isomorphic to the simplified target graph.

SI 首次定义检查方法：

> We reconstructed the building-unit connectivity graph from the assembly bond list and compared it with the target graph after simplifying both graphs.

随后可补具体判据：node count、edge count、degree sequence 和 graph isomorphism。审计其余位置 `topology recovery` 统一简写为 `building-unit connectivity check`。这一改法直接报告现有证据，保持期刊正文简洁；无需在摘要展开所有尚未实施的检验。

## 3. PLD/LCD：已经确认解析器列序颠倒

### 原始输出与解析器

Zeo++ 作者仓库文档说明 `-res` 的三个数依次为 maxIncDiam（最大内含球直径）、maxDiam（最大自由球直径）、incDiam（沿最大自由球路径的最大内含球直径）。即 Di、Df、Dif，常规孔指标对应 LCD=Di、PLD=Df。[作者文档](https://github.com/mharanczyk/zeoplusplus/blob/master/documentation/README)

`/Users/jing/PyCharmMiscProject/CO2_N2/src/co2_n2_screening/zeopp.py:37–42` 的 `parse_res_output` 当前返回 `{"pld": values[0], "lcd": values[1]}`，把前两列含义写反。正确映射应为 `{"pld": values[1], "lcd": values[0]}`。

本轮读取了全部 30 个候选的原始 `.res`：其中 9 个来自本仓库 targeted_zeopp 目录，其余 21 个来自相邻 CO2_N2 的 batch_002_clean_pool_5000 输出。逐项确认记录的 PLD 等于第一列、LCD 等于第二列；抽取实际解析函数并在原始文本上复现相同颠倒。故该问题已从疑点变为确认的解析错误。

例：`ffc_co2n2_cleanpool_01889.res` 三个数是 `27.57361 26.73585 27.57361`。正确 PLD=26.73585 Å，LCD=27.57361 Å。

### 影响和交付证据

- 本文 30 候选的两列数值及主表 PLD/LCD 范围应交换对应；组内跨度范围相应为 PLD 0.03–1.02 Å、LCD 0.18–0.92 Å。
- 原始 `.res` 无需重算。正式修复应改解析器后重新导出描述符及表格，保留更正记录；GCMC 原始吸附载量及从其计算的 S/WC/R/APS 不因两种孔描述符更名而数值改变。更广范围的筛选决策依赖另行核查，未在本轮推断。
- `match_pld_lcd_raw_audit_20261110.csv`：全部 30 候选的旧字段、正确字段和原始文件路径。
- `match_pld_lcd_corrected_group_preview_20261110.csv`：九组更正范围预览，尚未覆盖正式数据或正文。

## 4. 模板表与筛选节点库：两个清单的差异

正文 `tab:building_units` 列了五种 node templates 和一种 linker template；其标题声称覆盖 benchmark、construction、screening calculations。

30 个吸附候选的源表 `tables/ffc_co2n2_ordering_check.csv` 中，composition 字段实际共涉及十种节点类型。

| 正文表已列的五种节点 | 吸附 composition 另出现的五种节点 |
|---|---|
| 3c_B_Ch | 3c_pyrrole_Ch |
| 3c_triazine_Ch | 4c_In_1_Ch |
| 4c_B_Ch | 4c_Ni_1_Ch |
| 4c_Cd_1_Ch | 4c_V_1_Ch |
| 4c_Cu_1_Ch | 4c_tetrathiafulvalene_Ch |

最直观的例子是 G9 的三个模型：组成包含 `3c_pyrrole_Ch:4,4c_In_1_Ch:1,4c_Ni_1_Ch:1,4c_tetrathiafulvalene_Ch:1`，四种节点类型都没有列入该正文模板表。

这里发现的是论文模板清单与其声明覆盖范围不匹配；没有因此判定模拟使用了错误节点或组成。建议给现有表补五行并标注 CO2/N2 adsorption screening 用途；若表只打算列基准库，应收窄标题，并在 SI 单列吸附节点库。

## 核查状态

主稿 SHA-256 保持 `05270b0be8a26f3200a339a1d7ed0a43117f6a2bca17e9bcc5ef9de31b9a5102`。本轮未进行 LaTeX 编译或科学模拟；对 30 个原始输出的解析映射、字段一致性和组范围进行了核查。
