# MATCH 第二轮术语审查：剩余项目与代码核对

日期：2026-11-10。主稿校验和：`452ba52173c83c6e9907761b4d59b1c61848dd8ecb38ce856a420811ae56076e`。本轮审查保留已修正的 PLD/LCD、附录十节点模板表，以及作者明确要求保留的 topology recovery；不修改主稿。

执行更新（2026-11-11）：作者已批准并落实 R2-02、R2-03，涉及 Results 指派空间、覆盖核查、CASM 表注及结论轨道计数，中文同步。其余建议尚未获本轮修改指令；摘要词频及概念对应另见 `MATCH_ABSTRACT_TERM_AUDIT_2026-11-11.md`。

沿用 nature-polishing 技能，按算法型论文、英文全文及 MATCH 语境判断。本文在术语校准阶段；以下建议不等于已完成逐句语言终审。

## 一、核对后保留的表述

| 对象 | 核对依据 | 决定 |
|---|---|---|
| node assignment / coloring | 2.1 已将节点标签映射与其图着色表示对应；$X$ 独立按度数池选标签 | 保留，不机械统一成单一词，也不恢复作者已回滚的着色解释段 |
| ordering / arrangement | 讨论节点类型的周期空间排布，未被定义为热力学有序相 | 物理排布叙述中保留；算法输出按实际对象写 representative |
| fixed-point domain | $\mathrm{Fix}(H)$ 是被 $H$ 保持的指派集合，并用作搜索域 | 语境适当；不是数值固定点迭代，无需因 domain/set 两词出现而判错 |
| subgroup lattice / descent | 2.2 已明确子群包含偏序，且相邻任务可以不可比 | 方法名和解释可保留，不需要把 lattice 改成其他词 |
| Kabsch、$p_i$、$a_i$ | `compute_scaled_node_and_edge_vectors` 将邻接顶点位置放入 nvecs；`place_and_calculate_nodes` 用 `nvec-cvec` 和模板连接位点质心相对向量，调用 `minimize_rms_superimpose(...,center=False)`；核心 SVD 作行列式修正 | 与正文旋转最小二乘公式一致，保留。行/列向量写法对应转置约定，不是概念错误 |
| 全双射及最小 RMSD | 工具函数对不超过七个位点遍历全排列，严格小于当前最优值才更新；本文十种节点均为三/四连接 | 全双射和固定顺序破同分表述在本文节点库内有依据。不要扩展成任意高连接模板也遍历全部位点 |
| stabilizer fraction | 公式、匹配分组句和组成内比较的限定相容 | 保留；symmetry index 可精化但不是计算错误 |
| RMSD / closure error | 代码对连接位点实际间距与名义桥长的偏差计算均方根；SI 已定义 closure error | 名称可保留，不必另创替代指标名；重点是区别两种 0.8 Å 判据 |
| working-capacity fraction $R$ | $R=WC/q_{ads}$，正文已有明确定义 | 保留，不因原始字段叫 regenerability 就强制改名 |
| APS 的均值定义 | 表和 SI 区分 mean S × mean WC 与每次乘积的标准差 | 保留，不将二者合并为一个未注明的平均规则 |
| proper operations | 晶体学语境可指保向操作，已有 $\det R=1$ 限定 | 不将这个词直接判为错误；若改为 orientation-preserving，只属于范围更明确的措辞选择 |
| $\mathrm{Fix}(H)\cong\prod_P\mathcal C_{\deg(P)}$ | 此处为集合的自然一一对应，未赋予群结构 | 可保留现式；不应把原稿解读为已经错误声称群同构。必要时用文字说明 bijection 即可 |

## 二、建议优先修改的局部表述

### R2-01：连接标签 X 与占位元素 X 应分开

位置：主稿 444–463 行，SI 构建协议表的 Cleanup and export 行。

原句：

> Residual auxiliary \texttt{X} sites are then removed.

核对结果：

- 十种节点的连接位点标签都采用 X，但这些位点在 CIF 的 `_atom_site_type_symbol` 列有实际元素类型。九种模板的 X 位点为 C；`3c_pyrrole_Ch` 的 X1 为 N、X2/X5 为 C。
- 例如 `nodes/3c_triazine_Ch.cif` 的 X1、X3、X5 均是 C，`nodes/3c_pyrrole_Ch.cif` 明确写出 `X1 N`。
- `periodic_structure_analysis.py:519–548` 按真实元素和编号生成原子名称，同时保留原始连接标签用于追踪。
- `assembler.py:745–757` 按原子名称恢复 species，随后执行 `final_struct.remove_species(["X"])`。该删除条件是元素种类 X；带 X 连接标签的 C/N 不会因此被全部删除。端点融合时移除的匹配节点原子是另一操作。

建议在首次引入 X 时写：

> The label \texttt{X} marks connection sites independently of their atomic species.

清理句改为：

> Any remaining dummy atoms of species \texttt{X} are then removed.

中文分别为“标签 X 标记连接位点，与其原子元素种类区分”和“随后移除元素种类为 X 的残余占位原子”。协议表的 `Remove residual X` 同步写为 `Remove residual dummy X atoms`。这属于澄清代码已有语义，不改变清理规则或模板成分。

### R2-02：商图本身与其上的指派空间

位置：467 行。

原句：

> Their quotient-graph spaces span small exhaustive cases and larger cases tested with bounded enumeration.

判断：当前数学对象是商图 $Q$、指派空间 $X$ 和商空间 $X/\Gamma$。quotient-graph spaces 将图与其上定义的空间压成了一个不明确的名词组合；所指是不同规模的指派空间。

建议：

> The assignment spaces on their quotient graphs include small cases evaluated exhaustively and larger cases tested with bounded enumeration.

中文：“这些商图上的指派空间包括采用穷尽枚举的小规模算例，以及采用受限枚举检验的较大算例。”

### R2-03：轨道计数与原始指派计数

位置：CASM 表注 487 行、结论 663 行；472 行的覆盖核查句可一起校准。

原句分别称 `Orbits (Γ) are node-assignment counts` 和 `assignment counts depend on the acting group`。

判断：固定相同 $X$ 时，原始指派数 $|X|$ 不因作用群改变；CASM 比较改变的是等价类/轨道数 $|X/\Gamma|$。上下文能推断含义，但核心结论应直接说出计数对象。

建议表注：

> The Orbits ($\Gamma$) column reports the number of node-assignment orbits under the orientation-preserving permutation group used here.

建议结论：

> For \texttt{fjh}, the CASM comparison shows that the number of assignment orbits depends on the acting group.

472 行现称“每个指派由恰好一个候选轨道表示”，可以改为直接报告算法覆盖检验：

> For these finite spaces, we checked that the returned representatives covered $X$ without duplicate $\Gamma$-orbits.

这些改动均不改变当前计数、群或完整覆盖条件。

### R2-04：ordering energies 在此处可改得更具体

位置：623 行。

原句末尾：`to compare ordering energies at fixed composition`。

判断：这里是固定坐标、固定节点组成的单点模型能量比较。ordering energies 可以理解为“不同排布的能量”，并非必然用错；但也可能使读者联想到相对于参照态定义的有序化能。本段实际比较对象可以直接说清。

建议完整句：

> A separate set of constructed \texttt{ffc} and \texttt{fjh} models underwent single-point PBE+D3 calculations to compare their energies at fixed node composition.

随后保留原有“固定原子坐标”和 $\Delta E$ 的定义，不扩写稳定性结论。

### R2-05：secondary building units 对本文节点库的覆盖不够明确

位置：引言 139 行：`for its secondary building units and linkers`。

模板库中同时存在 triazine/pyrrole 等有机节点和金属节点。SBU 的具体化学定义取决于采用的分解方式；不能仅凭该词出现就断言所有有机节点都不允许叫 SBU。本文已经定义了更直接的 node 类型和节点模板，因此这里用 node building units 能无歧义地覆盖整个输入库。

建议短语：

> for its node building units and linkers

中文对应“节点构建单元与连接体”。本项是适用对象的精化，不是要求全稿删除 SBU 术语。SBU 化学含义的核对背景见 [Tranchemontagne 等，Secondary building units, nets and bonding](https://pubmed.ncbi.nlm.nih.gov/19384437/)；该文讨论特定金属羧酸簇，不能拿它推导所有文献对 SBU 的用法都只有一种。

### R2-06：proper pairwise joins 的修饰关系

位置：787 行。

原句：`Quota setup starts from cyclic subgroups and repeatedly takes proper pairwise joins.`

代码 `enumerator_rho_quota.py:90–99` 对生成元作并生成，并在 `combined.order() < group_order` 时加入。限制对象是并生成后仍为 $\Gamma$ 真子群的结果。

建议：

> Quota setup starts from cyclic subgroups and repeatedly forms pairwise joins that remain proper subgroups of $\Gamma$.

这保留“并生成”这一准确数学操作，不换成普通集合并集。

## 三、适合统一落地的较小项目

| 位置 | 当前形式 | 判断及建议 |
|---|---|---|
| 825–853 行 | $s$ 同时为稳定子阶数和晶胞扩张因子 | 两节有各自定义，尚未导致公式错误；建议扩张因子局部改为 $\lambda$，同步该表的全部阈值表达式，保留枚举 $s$ |
| 288 行，之后 DFS 表头 | depth-first search 未给缩写 | 首次改为 depth-first search (DFS)，数学内容不变 |
| 626 行 | symmetry index | 若统一到首次定义的量名，可写 stabilizer fraction $\rho$；当前不是另一个未定义计算指标 |
| 816 行 | stabilizer size | size 可表示基数，数学上成立；与主体一致可写 stabilizer order |
| 529–532、799 行附近 | layer 与 stratum | 搜索任务的阶数阶段与候选层分别指明；不因同义表达出现就指控概念错误 |
| 656 行图注 | low/medium/high strata | 区间分组本身也是一种分层；写 $\rho$ ranges/bins 能与精确 $\Omega_s$ 层区分，属于可选精化 |
| 877 行 | 500/60 Ry cutoffs | 数值正确，`dft_validation/config.py:103–104` 为 500/60，runner 写入 MGRID 的 CUTOFF/REL_CUTOFF。建议展开为 `a plane-wave cutoff of 500 Ry and a relative cutoff of 60 Ry` |
| 877 行 | SCF settings of diagonalization | 首次展开 self-consistent field (SCF)，并将搭配改为 `SCF calculations used diagonalization ...`；随后逐句润色时拆开当前长句 |

CP2K 的两个 cutoff 含义不同：CUTOFF 对应最细网格的平面波截断，REL_CUTOFF 控制 Gaussian products 的网格分配。数值不需要修改。[CP2K MGRID 官方文档](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/MGRID.html)

## 四、接续状态

- 本轮没有修改主稿、代码或结果表；不需要重新编译或运行科学计算。
- 已核对剩余术语所在的数学方法、构建实现、结果及 SI 句子，并确认 Kabsch 向量和本文三/四连接节点的全排列范围。
- R2-01、R2-02、R2-03 优先于纯措辞变化；R2-04 至 R2-06 可局部处理。其余项目应与后续逐句润色合并，避免为统一而改变本来准确的概念。
- 后续逐句稿应为每句记录“保留、可选精化、需要修改”及原因。本轮没有把主稿标记为语言终审完成。
