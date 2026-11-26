# 摘要与后文的术语一致性及单次用词审查

日期：2026-11-11。范围：`paper/main_MATCH.tex`。本轮按作者批准修改第二轮报告的 R2-02、R2-03（指派空间与轨道计数），摘要保持原文。本报告回答摘要用词是否与后文对应、是否仅用一次；不把原样短语出现次数直接等同于概念是否得到解释。

作者后续决定：摘要和贡献项均保留 `subgroup lattice`，已将两处表述落实为 `organize fixed-point searches through the subgroup lattice` 并同步中文。下列次数和“摘要未改”均描述此前审查快照；关于删去 lattice 的建议由作者本次决定取代。

## 1. 统计范围与方法

- 脚本：`audit_match_abstract_terms_20261111.py`。
- 数值证据：`match_abstract_term_counts_20261111.csv`、同名 JSON。
- 后文拆分为主文和 SI；包括英文正文、标题、图表说明和算法，排除中文/其他注释；不载入编译后参考文献条目。标题指正文/图表标题，不包含摘要前的论文题名。
- 统计不区分大小写，要求词或短语边界；去除 LaTeX 命令并将 CO$_2$/N$_2$ 归一化为 CO2/N2。单复数、词形和同义表达不自动合并。
- 表中“后文”不含摘要自身。短语之间有包含关系，计数不可相加推导总词数；后文还须人工检查概念对应。

## 2. 核心用词的对应情况

| 摘要用词 | 摘要次数 | 后文原样次数（主文 / SI） | 判断 |
|---|---:|---:|---|
| building-unit arrangements | 1 | 1 / 0 | 与引言首段一致 |
| node types | 1 | 1 / 1 | 另有单数 node type 及数学标签定义，所指一致 |
| fixed MOF net | 1 | 2 / 0 | 固定网络对象一致；详细方法还用 embedded net |
| periodic models | 2 | 8 / 1 | 反复使用，不是摘要孤立命名 |
| degree-compatible node assignments | 1 | 2 / 0 | 与引言、结果一致；2.1 给出度数池约束 |
| fixed linker template | 1 | 1 / 1 | 与构建输入及附录一致；固定在每条边上，未扩展成连接体枚举 |
| orientation-preserving | 1 | 6 / 0 | 与所选晶体学操作和置换群定义一致 |
| candidate equivalence | 1 | 2 / 0 | 与方法节名称及定义一致 |
| canonicalization | 1 | 1 / 0 | 与 2.2 叶节点全群规范化一致 |
| orbit representatives | 1 | 1 / 0 | 后文还有 canonical representative(s)，由其定义明确关联 |
| Kabsch-based node alignment | 1 | 1 / 0 | 与引言第二项贡献完全一致；2.3 有数学定义 |
| node alignment | 1 | 2 / 1 | 与流程图及构建协议一致 |
| linker placement | 1 | 2 / 0 | 与贡献项及流程图一致；2.3 展开操作 |
| oriented periodic models | 1 | 4 / 0 | 主文多处使用，方法说明了取向处理 |
| selected representatives | 1 | 4 / 0 | 与选择代表进行模型构建的限定一致 |
| geometric construction audit | 1 | 1 / 0 | 对应 Results 小节标题；未另创验证任务 |
| target topology | 1 | 1 / 0 | 对应审计结果；按作者决定保留此措辞 |
| short-contact check | 1 | 2 / 0 | 与 Results 审计一致；SI 给出接触判据 |
| reloadable | 1 | 1 / 0 | 与 Results 的 sanitized, reloadable CIF 一致 |
| CO2/N2 | 1 | 4 / 8 | 与吸附结果及计算协议一致 |
| adsorption screening | 1 | 1 / 6 | 与结果小节及模板用途一致 |
| composition groups | 1 | 1 / 0 | 对应 DFT 图注 eight composition groups；正文还有 fixed composition / within each composition |

这里的大多数术语在摘要内只出现一次，但后文也使用和解释。摘要无须为每个术语重复定义；也不应为减少重复引入另一套术语。

## 3. 原样仅在摘要出现一次的表达

| 原样短语（全文一次） | 后文对应表达/对象 | 判断与建议 |
|---|---|---|
| tunable structures | 引言 tunable porosity、diverse chemical functionality | 两者范围不完全相同。摘要以较宽的结构可调性概括背景，可保留；若换成 tunable porosity 是口径收窄，不能仅因前者一次出现而替换 |
| ordered arrangements | 引言 ordered patterns；后文 orderings / spatial arrangements | 所指均为节点类型空间排布，未发现概念冲突。可保留；若希望固定词组，可把引言 ordered patterns 对齐 ordered arrangements，属于可选统一 |
| selected symmetries | 引言 chosen symmetry group；下一句段落明确 orientation-preserving net symmetries | 首段问题句先概括、次段给限定，逻辑成立。可保留，不必每次重复全群定义 |
| compatible building-unit connections and orientations | 引言 compatible connection geometry and orientation；2.3 节点配准/连接体放置 | 复合句组非固定专名，后文有对应操作。连接与取向在语境中相关且各有所指，可保留 |
| quotient-graph colorings | 引言 quotient-graph coloring problem（单数）；2.1 graph colorings 与 $Q$、$X$ 的定义 | 单复数和句法形式不同，同一数学概念。保留复数 colorings 合适，不必为精确匹配频次而改 |
| subgroup-lattice fixed-point searches | 引言 deterministic subgroup fixed-point searches；2.2 Enumeration by Subgroup-Lattice Descent 和固定点域搜索 | 语义有支撑，但三重复合修饰比后文实际名称更重。建议摘要用 subgroup fixed-point searches，与贡献项一致；subgroup lattice 的算法组织放在方法节说明 |
| without duplication | 2.2 每个规范代表只记录一次；结果无重复轨道核查 | 属结果性质的普通表述，不是独立术语，可保留。duplication 与 duplicate 的词形不同不构成命名冲突 |
| crystallographic information files | 后文 CIF、CIFs | 全称仅在首次出现处展开正常；CIFs 在后文原样出现一次，单数 CIF 还在多处使用。保留全称和缩写 |
| single-point density functional theory | 后文 single-point DFT、PBE+D3 和计算协议 | 首次展开全称、后文使用缩写正常。DFT 本身在主文出现 4 次、SI 2 次（不含摘要），不存在只提一次的概念问题 |
| energy comparisons | 后文 comparison、relative energy、$\Delta E$ 和 energetic comparison | 同一动作及指标得到明确展开；不是另一个未定义能量。可保留。上一轮关于 ordering energies 的精化建议仍是单独事项，未经作者批准不修改 |

本表中的“全文一次”仅指精确短语，不表示其各组成词都只出现一次。当前最值得主动统一的是 subgroup-lattice fixed-point searches；其余项目大多属于语法变体、全称与缩写、或正常的背景概括。

## 4. 九句摘要逐句核对

| 句号 | 句子切入点 | 术语核对结果 |
|---|---|---|
| A01 | Metal-organic frameworks ... | MOF、应用范围和 building-unit arrangements 均有后文对应；tunable structures 较概括但没有概念错位 |
| A02 | Assigning node types ... | node types、fixed MOF net 准确；ordered arrangements 对应后文 orderings，未引入额外有序相定义 |
| A03 | For such a net ... | assignments inequivalent、selected symmetries、periodic models 与研究任务一致；连接和取向由构建步骤展开 |
| A04 | We formulate ... | degree-compatible、quotient-graph colorings、fixed linker template 符合 $X$ 定义和固定输入 |
| A05 | Selected orientation-preserving ... | 与所选保向操作诱导的 $\Gamma$ 定义一致；没有声称采用全部图自同构或原子模型空间群 |
| A06 | Subgroup-lattice fixed-point searches ... | 概念成立；建议缩短为 subgroup fixed-point searches 对齐贡献项，canonicalization 与不重复输出表述保留 |
| A07 | Kabsch-based node alignment ... | 与贡献句逐项相同；Kabsch 修饰节点配准，连接体放置另列，不存在作用范围挪用 |
| A08 | In the geometric construction audit ... | audit、target topology、short-contact、reloadable 与结果句对应；CIF 全称一次展开合理。遵从作者决定保留 topology recovery 相关措辞 |
| A09 | We also used models ... | CO2/N2、adsorption screening、single-point DFT、组成内比较有对应结果；未把几何审计与后续计算样本写成同一集合 |

## 5. 全文仅一次的单词形式

按脚本的英文词形统计，下列词形在摘要各一次，后文未原样出现：

`offer, motivating, systematic, produces, seek, orientations, duplication, density, functional, theory`

- 前五项是普通动词/修饰语，未承担新技术概念，不需要为了重复而替换或在后文补写。
- orientations 的后文形式为 orientation/oriented；duplication 有 duplicate 及去重机制对应。这是词形变化。
- density、functional、theory 出现在 DFT 全称里；后文采用 DFT 缩写。这是正常首次展开。
- 该清单不包括仅一次的全部长词组，也不是机械删词清单。词频工具不能替代上述概念和语境核对。

## 6. 本轮落实与边界

- R2-02 已落实：Results 的 quotient-graph spaces 改成 assignment spaces on their quotient graphs，并同步中文注释。
- R2-03 已落实：覆盖核查句明确 returned representatives 覆盖 $X$ 且无重复轨道；CASM 表注明确 node-assignment orbit counts；结论明确 number of assignment orbits，均同步中文。
- 摘要九句与本轮修改前逐字一致；本报告的摘要建议尚未写入主稿。X 清理句、SBU、ordering energies、proper pairwise joins 等未获本轮批准的建议仍未修改。
