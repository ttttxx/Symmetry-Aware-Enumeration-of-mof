# 逐项讨论：算法与论文一致性（2026-11-23）

2026-11-24更新：作者已批准并落实3/10/11/12/19/20/23/25/26/27/28/29/30/31，英文及中文同步。6/13/14/15/22/24明确不改，7/9此前也已决定不改。下文保存原讨论过程，当前状态以本段、主稿和MATCH_APPROVED_REVIEW_EDITS_2026-11-24.json为准；行号为历史快照位置。

阅读建议：先看6/14/15/22/24的简短解释，再一次看一个“改前—改后”条目。无需一次处理完。

| 项目 | 状态 |
|---|---|
| 1 | 执行：统一新标签池并重跑 |
| 2 | 执行：主表和附录使用当前实现新时间 |
| 3 | 2026-11-24批准并已应用 |
| 4 | 执行：配额及匹配前缀重跑 |
| 5 | 执行：正文引用有序标签池 |
| 6 | 2026-11-24明确不改 |
| 7 | 不改 |
| 8 | 执行：CLI/API默认auto，优先C++ |
| 9 | 不改 |
| 10 | 2026-11-24批准并已应用 |
| 11 | 2026-11-24批准并已应用 |
| 12 | 2026-11-24批准并已应用 |
| 13 | 2026-11-24明确不改 |
| 14 | 2026-11-24明确不改 |
| 15 | 2026-11-24明确不改 |
| 16 | 执行：补G1–G9组成映射 |
| 17 | 执行：用稳定子阶数替换symmetry index |
| 18 | 执行：补新基准包、图表映射和manifest；不推送远端 |
| 19 | 2026-11-24批准并已应用 |
| 20 | 2026-11-24批准并已应用 |
| 21 | 执行：去掉模糊This及推论连接 |
| 22 | 2026-11-24明确不改 |
| 23 | 2026-11-24批准并已应用 |
| 24 | 2026-11-24明确不改 |
| 25 | 2026-11-24批准并已应用 |
| 26 | 2026-11-24批准并已应用 |
| 27 | 2026-11-24批准并已应用 |
| 28 | 2026-11-24批准并已应用 |
| 29 | 2026-11-24批准并已应用 |
| 30 | 2026-11-24批准并已应用 |
| 31 | 2026-11-24批准并已应用 |

## 几个先说明白的问题

- **为什么还有旧计时？** 上一轮批准范围只覆盖方法，性能表尚待确认。现已明确统一新池、新实现，重跑后替换论文中的旧计时；历史原始文件仅作为过程档案保留。
- **6，不改的影响：** 不影响数学证明，也不影响所有层均为20,000的实验。差别是论文允许每层填不同数，代码只接受一个通用数。作者于2026-11-24决定不改，保留现有统一配额接口及论文表述。
- **14，说人话：** 论文要求输入对称操作彼此组合后仍在集合内，并且不会把三连接位点换到四连接位点。三个实验网络已核对通过。换新网络时，程序没有自动核对这些条件。现有实验不受影响；作者于2026-11-24决定不改通用入口。
- **15，不改的影响：** 已有30个模型的数值和组内比较不变。但别人不知道最初筛了哪些候选、如何挑选，难以复现样本来源，也不能判断30例能否代表更大空间。当前论文没有宣称随机代表性，所以主要是来源透明度问题。
- **22，不改的影响：** 句子位于引言，读成设计要求是合理的。风险是must被理解为每个指派一定能建出来。它不影响公式和实测结果；可保留，若改只需明确“目标要求”，不增加防御性限制句。
- **24，可以不改：** 粗分类表已覆盖正文提到的方法，维持现状。

## 用词约定

- Γ：节点指派判等的保向置换群；F：CASM原生完整因子群。
- enumeration wall time：包含任务准备、搜索、结果保存、解码与排序；native-core time只是其中一部分。
- stabilizer order / ρ：稳定子阶数 / 稳定子比例。精确层与低中高区间分开。
- 连接标签X与虚拟元素X：候选11项明确区分；未批准前不改主稿。

## 3. Results 3.1，原492–497行

**想解决什么：** 三类证据各说各的话：完整记录一致支持正确性；计时支持效率；F与Γ的计数差说明群选择的影响。不会写“结果一样，所以计时公平”，也不把小例较慢归因于未经拆分测量的单一成本。

**改前：**

```text
The \texttt{ffc} and \texttt{fjh} cases permit exhaustive enumeration over their degree-compatible node pools. For these finite spaces, we checked that the returned representatives covered $X$ without duplicate $\Gamma$-orbits.

For every exhaustive case, the represented orbit sizes summed to $|X|$. Distinct canonical representatives define disjoint $\Gamma$-orbits, so this sum counts the assignments covered and equality with $|X|$ confirms complete coverage. CASM provided a separate comparison on the same quotient graphs and node-label pools (Table~\ref{tab:casm}).
Its native full factor group $F$ includes orientation-reversing operations in addition to the orientation-preserving operations retained here. The \texttt{ffc} counts agree, providing a consistency check on the count computed for that benchmark space; this agreement does not establish that $F$ and $\Gamma$ define the same equivalence relation. For \texttt{fjh}, CASM gives fewer orbits under $F$, illustrating how the additional equivalences reduce the number of distinct candidates. The exact-coverage check therefore concerns the specified $\Gamma$ action, while the CASM comparison also exposes the effect of choosing a different equivalence group.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
For each exhaustive case, distinct canonical representatives defined disjoint $\Gamma$-orbits whose sizes summed to $|X|$, confirming complete coverage. The matched CASM workflow used the same assignment space, group action, and output fields. After normalization to a common canonical convention, the complete representative records agreed.

Enumeration wall times were lower for this method in 5 of the eight ffc and fjh cases in Table~\ref{tab:casm}. In the largest fjh case, 16,777,216 assignments yielded 1,078,816 orbits in 2.74 s, compared with 20.88 s for CASM. CASM was faster in the three smallest assignment spaces. These measurements support generating complete candidate libraries over the tested spaces, with relative performance dependent on the case.

CASM's native full factor group $F$ additionally includes orientation-reversing operations. Its orbit counts are reported separately in Supplementary Table~\ref{tab:orbit_group_counts}, showing the effect of the equivalence group.
```

**中文意思：** 先说明结果正确，再报告同群、同输出条件下的速度，最后把F群计数放回“群选择影响”的位置。数值已填入本轮重跑的实测结果。

**原拟稿的联动说明：** 主表与表注已获1/2项授权更新；本项新增解释段仍待确认。结论首段已有同群基准句，确认本项后补足其正文证据。

## 10. Algorithm 1，原357–359、363、386行

**想解决什么：** 让s和t的角色清楚，同时消除从k=0起步却未定义P₀的问题。搜索、合并和接收顺序都不变。

**改前：**

```text
\State build raw partition tasks $Q_0$ for the full group, identity and proper subgroups
\State order each partition's blocks by their smallest vertex index
\State form $Q$ from $Q_0$, merging proper-subgroup partitions; prioritize fewer vertex orbits, retaining construction order for ties
\For{each task $(\mathcal P,s)$ in $Q$}
\For{each admissible label for $P_k$ in the supplied order}
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
Build raw tasks $Q_0=\{(\mathcal P_H,|H|)\}$ for the full group, identity, and proper subgroups.
For each task, write $\mathcal P=(P_0,\ldots,P_{m-1})$, ordering blocks by their smallest vertex index.
[Keep the existing task merging and ordering steps.]
For each task $(\mathcal P,s)$ in $Q$, with $s$ the source subgroup order.
For each admissible label for $P_k$ in the supplied order.
```

**中文意思：** 每个任务保存“顶点分组＋来源子群阶数”。块从P₀开始编号。叶节点算出的实际稳定子阶数继续记为t。

**原拟稿的联动说明：** Algorithm 1的任务创建行和块顺序行；不改公式，不重写其他算法步骤。

## 11. 2.3，原464、482行

**想解决什么：** 标签是名字，元素是原子种类。代码按元素删除，不会因为C或N的连接标签叫X就把它删掉。

**改前：**

```text
Node and linker templates retain their specified conformations, with the auxiliary label \texttt{X} marking intended connection sites.
Residual auxiliary \texttt{X} sites are then removed.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
Node and linker templates retain their specified conformations. Connection-site labels beginning with \texttt{X} mark the atoms used for attachment; they retain their template chemical species during placement.

Residual dummy sites with chemical species \texttt{X} are then removed.
```

**中文意思：** X开头的标签标记连接原子；这些原子仍保留自己的元素。最后删除的是元素为X的虚拟位点。

**原拟稿的联动说明：** 2.3首次介绍X及清理句；操作表的Cleanup行相应写成“Remove residual dummy species X”。其余X–X连接记号保留。

## 12. 2.3，原479行

**想解决什么：** 顺序直接按assembler.py执行顺序写，明确审计量取自融合前。这里只改说明，不改数据。

**改前：**

```text
After these rotations, node-side and linker-side \texttt{X} endpoints with a Cartesian separation strictly below $0.8\,\text{\AA}$ are merged one-to-one. After fusion, the remaining topology-linked node--linker \texttt{X}--\texttt{X} pairs define the bridge connections. Their CIF bond-table entries use a nominal length of $1.485\,\text{\AA}$ without displacing the atomic coordinates; the geometric audit measures the actual \texttt{X}--\texttt{X} separation before this nominal value is applied. The numerical criteria used in these steps are listed in Supplementary Section~\ref{app:reconstruction_protocol}. The fixed-geometry scan in Appendix~\ref{app:x_fusion} examines this conservative operational cutoff.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
After these rotations, the geometric audit measures the actual separations of topology-linked node--linker \texttt{X}--\texttt{X} pairs before endpoint fusion. Node-side and linker-side \texttt{X} endpoints with a Cartesian separation strictly below $0.8\,\text{\AA}$ are then merged one-to-one. After fusion, the remaining topology-linked pairs define the bridge connections. Their CIF bond-table entries use a nominal length of $1.485\,\text{\AA}$ without displacing the atomic coordinates. The numerical criteria are listed in Supplementary Section~\ref{app:reconstruction_protocol}. Appendix~\ref{app:x_fusion} reports the fixed-geometry scan of the fusion cutoff.
```

**中文意思：** 轴向旋转后先测真实端点间距，再融合端点。剩余桥键在CIF键表中写名义长度，坐标不跟着移动。

**原拟稿的联动说明：** 2.3本段；附录Construction audit protocol可在closure error定义中补“before endpoint fusion”，与正文一致。

## 13. 3.2，原559行

**想解决什么：** RMSD和95.7%描述同一批数据，放在一起就不会让人把融合距离与闭合误差混淆。

**改前：**

```text
The node--linker bridge lengths quantify how closely rigid building-unit placement matches the target connection geometry. Across 3,888 topology-linked \texttt{X}--\texttt{X} pairs, the RMSD from the nominal $1.485$\,\AA{} bridge length was 0.388\,\AA{}. The fusion tolerance is 0.8\,\AA{}. Separately, the closure audit measures the actual \texttt{X}--\texttt{X} separation against the nominal $1.485$\,\AA{} bridge length; 95.7\
```

**历史候选，作者明确不改：**

```text
The node--linker bridge lengths quantify how closely rigid building-unit placement matches the target connection geometry. Across 3,888 topology-linked \texttt{X}--\texttt{X} pairs, the RMSD from the nominal $1.485$\,\AA{} bridge length was 0.388\,\AA{}. In the same closure audit, 95.7\% of separations lay within 0.8\,\AA{} of the nominal length. Together with topology recovery and the short-contact check, these measurements support using the constructed models as periodic starting structures for downstream calculations. Structural relaxation remains a subsequent step. The contact criterion and sample selection are detailed in Supplementary Section~\ref{app:reconstruction_audit}.
```

**中文意思：** 先说同一批3,888对端点的RMSD为0.388 Å，接着说95.7%的误差在0.8 Å内。融合阈值只留在方法和操作表。

**原拟稿的联动说明：** 只改3.2本段；方法479行和操作表的0.8 Å融合阈值保留。与12项的采样时点一致。

## 19. 引言第2段，原134行段首

**想解决什么：** 去掉“下面有几种策略”的目录句，内容和引用不增加；读者一开始就看到方法具体做什么。

**改前：**

```text
MOF candidate generation has followed several strategies. Topology-based assembly methods such as ToBaCCo place building units on selected nets~\cite{colon2017topologically}, while reverse-topological methods match unit geometry to topological constraints~\cite{wu2024precisionmof}.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
Topology-based assembly methods such as ToBaCCo place building units on selected nets~\cite{colon2017topologically}. Reverse-topological methods match unit geometry to topological constraints~\cite{wu2024precisionmof}.
```

**中文意思：** 直接从拓扑组装的具体操作开始，再说反向拓扑；后面的遗传算法和生成模型句子保留。

**原拟稿的联动说明：** 只改引言第2段开头；粗分类表维持原样。

## 20. 引言第4段，原139行首两句

**想解决什么：** 本段负责说明“为什么值得研究有序排布”。用同组成不同排布切入即可，度数相容和逐点指派定义留给方法。

**改前：**

```text
For a fixed MOF net, an analogous assignment specifies a degree-compatible node type at each site. Different assignments can produce ordered patterns, including patterns with the same composition.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
On a fixed MOF net, node types can form distinct spatial arrangements at the same composition.
```

**中文意思：** 固定MOF网络上，即使组成相同，节点类型也可以有不同的空间排布。

**原拟稿的联动说明：** 与已批准21项的明确主语句衔接；23项表后段负责研究范围，不再重复动机。

## 23. 引言表后段，原178行

**想解决什么：** 唯一职责：界定研究范围。前一段已说为何研究排布、为何需要判等和建模；比较表已说明本文位于哪类工作；后一段贡献项说具体得到什么。这里应只交代本研究固定哪些输入、变化哪个对象。因此删掉重复任务句，并把商图着色、Γ轨道的正式定义全部留给2.1。

**改前：**

```text
For a fixed embedded MOF net, the task is to identify symmetry-distinct node assignments and construct periodic models with compatible building-unit connections. Within a chosen repeat cell, we formulate degree-compatible node assignments as a quotient-graph coloring problem and hold one linker template fixed on each edge. Colorings related by selected orientation-preserving symmetries of the embedded net define one candidate orbit; selected representatives are then used to construct oriented periodic models.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
We study node-type arrangements within a chosen repeat cell of a fixed embedded MOF net. Linker templates are specified for the edges and held fixed throughout node-assignment enumeration.
```

**中文意思：** 我们研究固定嵌入MOF网络在所选重复晶胞内的节点排布。边上的连接体模板预先指定，在节点指派枚举中保持不变。

**原拟稿的联动说明：** 引言178整段与20/21项一起审视。若26项采用统一连接体范围，第2句同步改为“A single linker template is held fixed on every edge throughout node-assignment enumeration.”。贡献项和2.1本轮不据此重写。

## 25. 2.3段首，原420行

**想解决什么：** 直接进入本节新任务“Γ从哪里来”。不再先回顾Γ用于候选判等。后段接着定义这些操作和映射。

**改前：**

```text
The group $\Gamma$ defines candidate equivalence in the coloring problem. To apply this definition to a periodic framework, we derive its vertex permutations from the crystallographic symmetries of the embedded net.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
The crystallographic symmetries of the embedded net induce vertex permutations in the chosen repeat cell. We use their orientation-preserving permutation image as $\Gamma$.
```

**中文意思：** 嵌入网络的晶体学对称操作在所选重复晶胞中诱导顶点置换。我们把其中保向操作的置换像用作Γ。

**原拟稿的联动说明：** 只改2.3开头两句；紧接着的形式定义保留。

## 26. 输入范围；摘要120、引言178、方法199/464、流程图注205

**想解决什么：** 建议把范围收准到当前程序和实验：所有边都用同一种连接体。这样网络对称操作不会交换不同种类的边模板，也不需要新加一套群约束。这会收窄“每条边可分别固定不同模板”的一般化范围，但不改变现有实验。若作者确实需要不同边使用不同模板，则应另做边模板保持群的实现和论证，不能只改一句。

**改前：**

```text
We fix the embedded net, a chosen repeat cell, and one linker template for each graph edge.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
We fix the embedded net, a chosen repeat cell, and a single linker template used on every graph edge.
```

**中文意思：** 固定网络和重复晶胞，并在所有边上使用同一种固定连接体模板。

**原拟稿的联动说明：** 摘要改“with a single fixed linker template on every edge”；2.1首句、2.3模板句、图1图注及23项的范围句一起统一。X、Γ公式和全部当前数据保持不变。

## 27. 2.1规范代表定义244及2.2配额保证344

**想解决什么：** 这两种顺序允许不同。把它们分清，既符合代码，也不需要改变现有定义或证明。新的统一计时池恰好按字典序输入，但算法不依赖这种巧合。

**改前：**

```text
The implementation compares label strings in lexicographic order.
If no task times out and $M_{\mathrm{total}}\geq\sum_sB_s$, stratum $s$ returns the first $\min\{B_s,|\Omega_s|\}$ representatives in the fixed task, vertex-orbit, and label order. Reaching $B_s$ fills the quota; completing all tasks for that stratum before filling the quota establishes exhaustion. The coverage and sequence guarantees are proved in Supplementary Section~\ref{app:quota_proof}.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
Canonical comparisons use the lexicographic order of label strings. DFS visits the admissible labels in their supplied pool order, which determines the traversal sequence within each task.

If no task times out and $M_{\mathrm{total}}\geq\sum_sB_s$, stratum $s$ returns the first $\min\{B_s,|\Omega_s|\}$ representatives in this fixed traversal sequence.
```

**中文意思：** 判定代表用标签字典序；DFS先试哪个标签由输入列表决定。配额截取的是固定遍历产生的接收序列。

**原拟稿的联动说明：** 2.1比较顺序句；2.2配额保证句；附录L_s定义用同一“supplied pool order”。Algorithm 1原supplied order保留。

## 28. 附录lil段，原880行

**想解决什么：** 表里已经有指派数、轨道数、时间。正文保留新增网络及其证据范围，减少重复报数。

**改前：**

```text
An additional \texttt{lil} case used $|\mathcal C_3|=4$, $|\mathcal C_4|=3$, and $|\Gamma|=4$. Its 5,308,416 assignments yielded 1,328,832 orbits. The wall times were 1.952\,s for this work and 6.817\,s for the matched CASM workflow. This case was used for enumeration only.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
An additional \texttt{lil} case used $|\mathcal C_3|=4$, $|\mathcal C_4|=3$, and $|\Gamma|=4$. This case tested enumeration on a further net; no atomic models were constructed for it.
```

**中文意思：** 额外lil算例使用4种三连接标签和3种四连接标签，群阶数为4。它补充另一网络上的枚举检验，不构建原子模型。

**原拟稿的联动说明：** 只改lil说明段；其计时和计数仍按本轮重跑结果放在附录表。

## 29. DFT图注，原673行末句

**想解决什么：** 一个精确ρ值才是配额方法的一层；这里合并成三个区间，因此用ranges避免同词两用。

**改前：**

```text
Vertical guides indicate low ($\rho<0.25$), medium ($0.25\le\rho\le0.75$), and high ($\rho>0.75$) strata.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
Vertical guides indicate low ($\rho<0.25$), medium ($0.25\le\rho\le0.75$), and high ($\rho>0.75$) $\rho$ ranges.
```

**中文意思：** 竖线标示低、中、高三个ρ区间。

**原拟稿的联动说明：** 只改图注英文及中文；阈值、图件和精确层的名称不变。

## 30. 结论第2段首句，原683行

**想解决什么：** 让每类证据承担它实际支持的结论，避免把“枚举正确”直接写成“几何构建正确”的支持。其余孔结构差异、排名波动和弛豫边界保持。

**改前：**

```text
Exact count checks, construction audits, and matched screening comparisons support generating periodic starting models from symmetry-distinct assignments for downstream screening.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
The 72-model audit established successful construction from discrete assignments. Separate matched screening sets demonstrated the use of constructed models for downstream evaluation.
```

**中文意思：** 72例审计确认离散指派成功转为周期模型；另行选取的匹配筛选集合展示构建模型的后续评价用途。

**原拟稿的联动说明：** 结论第2段首句替换为这两句；轨道覆盖已在结论首段交代。后续孔结构差异、排名波动和弛豫边界保持。审计和筛选样本明确分开。

## 31. DFT协议，原1051行

**想解决什么：** 这里交代报告结果采用什么节点就够了。未展示结果的Cu探针和视觉用途会让人以为还有一组要解释的DFT证据。

**改前：**

```text
The reported single-point DFT sets used the closed-shell Cd node; Cu served as an optional open-shell probe and also appeared in visual examples.
```

**已批准并应用的候选（联动排版以主稿为准）：**

```text
The reported single-point DFT sets used the closed-shell Cd node.
```

**中文意思：** 报告的单点DFT集合使用闭壳层Cd节点。

**原拟稿的联动说明：** DFT协议删去Cu旁支；模板表Cu行删去optional open-shell probes，只保留实际报告的构建、吸附和视觉用途。

## 已落实部分的核对结果

- 主表及附录：108次独立进程全部穷尽，9例完整记录与同群CASM一致；当前ABI-4、新字典序标签池。主表8例中5例更快，最大fjh约2.740秒对20.88秒。
- vba：12次当前C++测量；配额总枚举时间中位数约0.2366/0.2849秒。代表总数和各层数量不变；新三标签池的DFS中间层组成向量数为43，配额对应层为54，已更新表中括号数。
- 默认CLI/API改为auto，优先C++，不可用时回退优化Python；legacy保留显式入口。相关189项测试通过。
- 补齐G1–G9组成、统一稳定子术语、澄清引言This句，并完成本地数据包与图表映射。未推送远端。
- 主稿编译34页，新表目视检查通过；5处既有overfull保留，没有新增overfull。

## 2026-11-24落实结果

十四项批准内容已写入主稿并同步中文；13项保留原闭合误差段落。算法只澄清来源子群阶数和从0开始的块编号，未改变实现。所有数学公式、数值结果表及测量源码/输入哈希不变。历史候选中的待定语句不再作为新一轮授权或待办。
