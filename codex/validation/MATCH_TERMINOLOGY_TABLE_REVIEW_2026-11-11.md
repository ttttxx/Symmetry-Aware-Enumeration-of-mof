# MATCH 全文术语与表注复核

最新执行：作者明确确认后，仅落实吸附表单位入表头、删除 caption 单位说明、压缩分组/范围/平均口径，以及删除 vba caption 重复的 ρ 定义和耗尽/配额判据。运行参数与吸附跨度保留在各自 caption，正文不变；无表下 Notes。下文原始搬移方案已被否定，不作为当前状态。

更新：作者随后否定第 3 节的两表分工调整，并澄清原请求为讨论，不代表授权直接修改。两张表及相邻正文已恢复至上一轮修改前；运行参数不应堆入正文，吸附跨度应在 caption 简要概括，单位入表头作为待讨论方案。以下第 3、4 节关于两表改动的记录仅用于历史追溯，不能作为当前已接受状态。

日期：2026-11-11。对象：`paper/main_MATCH.tex` 全文（含摘要、算法、图表及附录）。沿用前两轮术语核对；本轮落实作者点名的表注、连接体说明和标签调整，其余术语建议列出供逐项决定。

## 1. 当前术语表与保留判断

| 对象 | 建议统一形式／当前定义 | 本轮判断 |
|---|---|---|
| 构建单元总称 | building unit | 主体使用稳定；引言一次 building blocks 可统一，属于措辞选择 |
| 节点类型指派 | node assignment / node-type assignment | 离散变量，区别于几何 node placement |
| 图上的表示 | coloring | 2.1 已与 node assignment 建立对应，保留 |
| 空间排布 | ordering / arrangement | 排布或其对应模型的语境合理，保留 |
| 算法返回对象 | canonical representative | canonical ordering representatives、canonical assignments 可进一步统一 |
| 候选与轨道 | candidate orbit / Γ-orbit | 同一 Γ 轨道作为一个候选，保持 |
| 原始计数 | raw assignment count，\|X\| | CASM 表 Naive 表头可改为带 \|X\| 的明确名称 |
| 等价约简后的计数 | orbit count，\|X/Γ\| | 已修正，未再发现将其与原始计数混淆的句子 |
| 所选作用群 | orientation-preserving permutation group Γ | 区分晶体学操作及其置换像；保留 |
| 搜索域 | fixed-point domain / Fix(H) | 可与 fixed-point set 对应，不必机械替换 |
| 方法结构 | subgroup lattice / subgroup-lattice descent | 已说明任务可对应不可比子群；按作者决定保留 |
| 代表所属精确层 | exact-ρ stratum / Ω_s | 层定义为精确稳定子阶数；区别于图示区间分组 |
| 稳定子阶数与比例 | stabilizer order s；stabilizer fraction ρ | symmetry index、stabilizer size 是可统一的别称 |
| 分层输出限制 | per-stratum quota B_s | 方法仍用 budget，算法及结果使用 quota，建议统一 |
| 全局限制 | global output cap；task time limit/timeout | 与 B_s 分开命名，现有边界保留 |
| 构建结果 | constructed periodic model / pre-relaxation model | construction 为首次建模；reconstruction 的历史标识符可保留，无需批量改标签 |
| 几何核对 | geometric construction audit；topology recovery；short-contact check | topology recovery 按作者既定决定保留 |
| 配准指标 | root mean square deviation (RMSD) | 已展开；节点向量配准与连接间距审计须各自写清参照量 |
| 孔描述符 | pore-limiting diameter (PLD)；largest-cavity diameter (LCD) | 当前列值及跨度沿用修正结果，二者未互换 |
| 吸附指标 | selectivity S；working capacity WC；working-capacity fraction R；adsorbent performance score APS | 公式与正文一致；R 不强行改成原始字段 regenerability |
| 平均口径 | APS = mean S × mean WC | 表注和附录相容；运行间标准差采用各次 APS 的定义另行保留 |
| 能量比较 | single-point energies at fixed node composition | 当前 ordering energies 可精确化，不能扩展为弛豫稳定性结论 |
| 模板标识 | 1B_4H_Ch 等原始名称 | 保留实际标识符，固定连接体说明现并入 Results 首段 |

## 2. 仍需处理的术语与符号

下列项目本轮未自动写入主稿，延续此前逐项敲定的流程。位置按本轮修改后源码，后续可能随编辑移动。

| 位置 | 当前形式／问题 | 建议与重要性 |
|---|---|---|
| 2.3 连接标签及清理句，约 444、462 行；附录构建协议表 | X 同时表示连接标签和占位元素；Residual auxiliary X sites 容易被读为删除全部带 X 标签的 C/N 位点 | 优先澄清：连接标签独立于真实元素；清理删除 species X 的占位原子。代码与模板证据沿用 ROUND2 的 R2-01 |
| 引言第四段，约 139 行 | secondary building units and linkers | 可用 node building units and linkers 覆盖有机节点和金属节点；属于适用对象精化 |
| 2.2，约 288 行 | depth-first search，后文直接 DFS | 首次补 depth-first search (DFS)；定义补全 |
| 2.2，约 310、321 行 | budget B_s / quota B_s | 同一输出数量约束，建议统一 quota，中文同步为配额 |
| 审计及 DFT，约 536、623 行 | canonical ordering representatives；canonical assignments | 算法对象建议统一 canonical representatives；模型和排布语境保留 ordering |
| DFT，约 620、623 行 | ordering energies；symmetry index | 分别建议 energies at fixed node composition、stabilizer fraction ρ；前者是精化，后者是名称统一 |
| DFT 图注，约 653 行 | low/medium/high strata | 建议 ρ ranges/bins，与精确 Ω_s 层区分，保留现有区间端点 |
| 附录配额证明与审计，约 796、813 行 | layer coverage；stabilizer size | 可统一 stratum coverage、stabilizer order |
| 附录构建协议，约 822 行及表格 | 扩张因子 s 与稳定子阶数 s 重名 | 两节各有定义，未发现数值错误；建议局部将扩张因子改为 λ 并同步所有阈值公式 |
| 附录配额证明，约 784 行 | proper pairwise joins | proper 应明确修饰并生成后仍为 Γ 真子群的结果，沿用 R2-06 建议 |
| 附录 DFT，约 874 行 | 500/60 Ry cutoffs；SCF 未首次展开 | 可区分 plane-wave cutoff / relative cutoff，并展开 self-consistent field (SCF)；不改数值 |
| 2.3 保向操作解释 | a reflection would reverse template handedness | 涉及手性的化学解释可限定 for a chiral template；属于范围说明，不改所选作用群 |

上述多项属于可选精化；不应据此将全文现有概念一律判为错误。优先项仍为 X 标签与元素歧义，其次是缩写定义及同一量的名称统一。

## 3. 已落实的局部调整

- Results：将固定连接体名称并入开头的输入模板说明，删除 Exact and External Enumeration Checks 中单独成段的 All reported models...；附录模板表继续保留对应行。
- 表 vba_bounded：运行参数移到表前；正文已有相同输出数、精确层及耗尽/配额填满判据，表注不再重复。表注仅保留读表所需的匹配关系、括号内计数含义与任务耗时定义。Layer time 改为 Task time，明确这是固定点任务的耗时。数值不变。
- 表 ffc_co2n2_ordering：表头增设单位行，PLD/LCD 为 Å，WC/APS 为 mmol g⁻¹，S/R 的单位记为 1（无量纲）；表注保留分组、最小—最大范围和跨运行均值定义。九组共 30 个模型及三个组内跨度移至表后正文。数值不变。
- Geometric Construction Audit：保留 sec:reconstruction，删除 sec:validation。两个不同 label 原本指向同一小节，LaTeX 允许这种别名写法；当前主稿无任一标签的引用，删除不影响交叉引用。无法仅凭源码确定历史上为何留下别名。
- 英文及中文辅助注释同步；未修改算法、数据、科学结论或此前作者决定保留的表述。

## 4. 验证与覆盖

- 已逐节读取标题/摘要、引言、2.1–2.3、算法 1–2、枚举/几何/吸附/DFT 结果、结论及附录 A.1–A.8，并对关键变体做无注释正文检索。
- 静态检查：无重复 label、无缺失的正文 ref/pageref/eqref 目标（LastPage 由宏包生成）；两张表从 midrule 到 end{tabular} 的全部数据行与修改前逐字相同。
- 修改前后均完整 LaTeX 编译，均为 29 页；目视核对两版第 15、17 页及整稿缩略图，新表头和表注未溢出。
- 原有六项 overfull hbox、算法 1 的 Float too large（16.13359 pt）和 Hfootnote.1 警告仍在；修改前后警告类型和幅度一致。本轮不作整稿排版修复。
- 全文术语复核已完成；其余建议尚未全部应用，不能据此将稿件标记为术语或语言终审完成。
