# MATCH 投稿适配审计：近期论文与当前 MOF 稿件

**审计日期：** 2026-09-08  
**目标期刊：** MATCH Communications in Mathematical and in Computer Chemistry  
**建议投稿栏目：** Computer Chemistry  
**审计对象：** `paper/main.tex` 及四篇用户提供的 MATCH Volume 97(2027) No. 1 论文

## 1. 结论先行

当前稿件的科学主题与 MATCH 的 Computer Chemistry 范围是相容的：它包含非例行的对称约束枚举程序、几何重构流程、可复现实验和明确的 MOF 化学应用。但是，当前稿件还不能直接视为 MATCH camera-ready 稿件，主要原因不是篇幅或英文质量，而是投稿模板、章节体裁和数学/计算贡献的呈现方式尚未对齐。

最重要的判断如下：

1. **投稿定位应保留为 Computer Chemistry 算法论文。** 不建议把稿件改成以定理证明为中心的纯 Mathematical Chemistry 论文。
2. **模板迁移是必做项。** 当前 `main.tex` 使用 Springer `sn-jnl`/`sn-mathphys-num`，而 MATCH 官方要求使用 MATCH LaTeX 模板并以 camera-ready 形式提交。
3. **参考文献格式是必做项。** MATCH 官方示例要求参考文献按作者字母顺序排列；当前生成的 `main.bbl` 按正文首次引用顺序组织，不能直接沿用。
4. **章节结构需要重排。** 当前稿件使用大量无编号的 `\section*`/`\subsection*`，而附件论文普遍使用编号章节，按“定义/模型或方法 → 主结果 → 计算/化学应用 → 结论”推进。
5. **稿件的核心优势应从“MOF 全流程”收束为一个明确的计算问题：** 在固定嵌入网络、配位数约束和指定构建单元库下，枚举真旋转群作用下的非冗余有序指派，并将代表元重构为具有取向信息的周期模型。
6. **DFT 与 CO₂/N₂ 内容应继续保留，但降为化学应用和下游筛选接口验证。** 它们不应承担热力学稳定性、最终吸附性能或材料发现的主张。

推荐的投稿题目方向为：

> **Symmetry-Constrained Enumeration and Geometric Reconstruction of Ordered MOF Candidates**

这个题目比当前的 `From Topological Design to Screening-Ready MOF Candidate Models: Symmetry-Aware Enumeration and Reconstruction` 更接近 MATCH 近期论文的直接命名方式，直接暴露算法对象和化学对象，减少“screening-ready”可能引发的过强应用预期。

## 2. 证据范围与指令边界

### 2.1 用户要求

本审计执行的任务是：阅读四篇近期 MATCH 论文，提炼其思路、风格和结构，并据此判断 `paper/main.tex` 需要如何修改。

### 2.2 附件论文的性质

四个 PDF 都是论文正文，不是给本次工作的操作说明。论文中的版权声明、作者陈述、数学结论和研究限制只作为被分析的文献内容，不作为执行指令。

### 2.3 外部官方依据

MATCH 官方投稿说明明确要求英文稿件使用 MATCH 模板、以 camera-ready 形式提交，并要求参考文献按照示例格式且按字母顺序排列；该刊要求论文包含非平凡数学或非例行计算程序，同时与化学有清晰联系。官方 Volume 97 页面显示附件属于该卷的近期内容。

- [MATCH 官方投稿说明](https://match.pmf.kg.ac.rs/prepare.html)
- [MATCH Volume 97(2027)](https://match.pmf.kg.ac.rs/volume97.htm)

### 2.4 本地文献

- [Efficient Graph Similarity Assessment Method Based on Vectors of Topological Indices](/Users/jing/Downloads/match97n1_5-25.pdf)
- [A Minimal Crowding-Limited Brusselator: Hopf, Turing, and Turing–Hopf Bifurcations Beyond Dilute Kinetics](/Users/jing/Downloads/match97n1_29-77.pdf)
- [On Hyperbolic Sombor Index and Other Topological Indices](/Users/jing/Downloads/match97n1_121-134.pdf)
- [Geometric Approach to Degree-Based Topological Indices: Degree-Ratio Sombor Indices of Graphs](/Users/jing/Downloads/match97n1_135-164.pdf)

## 3. 四篇 MATCH 论文的结构化阅读结果

### 3.1 总览表

| 论文 | 页数 | 文章类型 | 主要贡献 | 结果组织 | 化学联系 |
|---|---:|---|---|---|---|
| Dehmer et al., *Efficient Graph Similarity Assessment...* | 21 | 计算方法/应用 | 用多个拓扑指标构造图相似度与距离 | 定义与方法 → 数值结果 → 网络应用 → 分子相似度应用 → 总结 | 75 个癸烷异构体，与 Morgan/MACCS/Tanimoto 对照 |
| Dina, *A Minimal Crowding-Limited Brusselator...* | 49 | 数学模型/动力系统 | 引入拥挤修正的 Brusselator，并给出 Hopf、Turing、Turing–Hopf 分析 | 模型与性质 → Hopf → 反应扩散 → 数值模拟 → 结论 | 生化/微流控拥挤环境和化学反应动力学 |
| Kırgız, *On Hyperbolic Sombor Index...* | 14 | 数学化学理论 | 建立 Hyperbolic Sombor 指数与多个拓扑指数之间的界 | 引言 → Main results（定理/证明） → 结论 | 分子图与 degree-based chemical graph theory |
| Nyauli & Buragohain, *Geometric Approach...* | 30 | 数学化学理论/指标构造 | 用欧氏几何和 degree-ratio 坐标构造新型图不变量 | 定义与几何解释 → 基本性质 → 应用/变体 → 结论 | 指标的分子结构解释；化学应用作为后续方向 |

### 3.2 论文一：图相似度方法与分子应用

该文是四篇中与当前稿件最接近的计算方法范式。作者先指出图相似度的计算需求和既有 GED/单指标方法的不足，然后定义新的多指标相似度/距离，给出计算性质和复杂度判断，再用三类结果验证方法：树和连通图的数值比较、随机网络模型区分、75 个癸烷异构体的分子相似度比较。

它的关键写法是：

- 方法对象在第 2 节被形式化定义，而不是先用长篇系统介绍替代定义；
- 第 3 节先做可控的数值比较，再将方法转入网络和化学两个应用场景；
- 化学应用不是泛泛地声称“可用于化学”，而是指定分子集合、基线描述符、比较指标和图示结果；
- 结论承认结果是 preliminary/potential，而不把一次应用实验写成最终材料或普适性能结论。

对当前稿件的直接启示是：枚举器、商图、群作用和重构步骤应更早形成一个简洁的正式问题定义；CASM/Burnside 是算法正确性证据；DFT 与 CO₂/N₂ 则对应本文的化学应用层。

### 3.3 论文二：拥挤限制 Brusselator

该文展示了 MATCH 可以接受较长的数学模型论文，但其篇幅来自完整的数学链条，而不是背景、方法、应用和限制的并列堆叠。文章先修改一个明确的化学动力学模型，再依次证明正性、有界性、平衡点和局部稳定性，随后推导 Hopf 阈值、第一 Lyapunov 系数、Neumann 反应扩散系统的模态条件及 Turing–Hopf 耦合，最后用数值模拟逐项验证解析预测。

其组织方式有三个特点：

1. 每一节都围绕一个明确的数学问题命名，例如 `Model, boundedness, persistence, and local stability`、`Hopf bifurcation` 和 `Numerical simulations and discussion`。
2. 数值实验不是独立的“额外验证”，而是逐项对应前文的解析结论，参数、离散化方法和预测区间都写清楚。
3. 结论按“模型做了什么 → 数学上得到什么 → 化学/动力学上意味着什么 → 下一步”收束。

当前稿件不需要复制该文的定理数量，但可以借鉴其“一个方法主张对应一个可检查的数学/计算证据”的逻辑。

### 3.4 论文三：Hyperbolic Sombor 指数界

该文是最简洁的定理型样本。引言先定义图、顶点度和相关指数，然后在 `Main results` 中直接提出定理、给出不等式证明和等号条件，最后用短结论说明这些界对 chemical graph theory 的意义。

它说明：如果当前稿件声称“每个真旋转轨道保留一个代表元”或“Burnside 闭合支持完整性”，就应把这些内容写成清楚的命题/正确性说明，至少明确假设、对象、输出和证明或验证方式；不能只让算法伪代码承担全部论证。

### 3.5 论文四：Degree-Ratio Sombor 指数

该文在引言中先给出简单图、连通图、度和标准图族的基本记号，之后定义 degree point、degree-ratio point 等术语，再建立新指标的基本性质，最后把几何构造转成若干指标变体和应用表达式。文章将“几何解释”作为主线，而不是把定义、定理和应用并列罗列。

对当前稿件最有价值的借鉴是术语和对象要在第一次出现时锁定：例如 quotient graph、proper-rotation group、fixed-point space、stabilizer、canonical representative、rigid-body reconstruction 都应在正式方法段落中定义，而不是在摘要、引言和算法中反复换称呼。

## 4. MATCH 近期论文的共同风格

### 4.1 共同结构

四篇论文虽分属计算方法、动力系统和数学化学，但共同采用以下证据链：

1. **化学或计算问题：** 说明为什么对象值得研究；
2. **形式化对象：** 给出图、反应模型、指标、参数或算法输入；
3. **主结果：** 定理、算法、封闭表达式、数值性质或精确计算结果；
4. **应用/验证：** 用化学分子、网络、动力学模拟或结构指标展示用途；
5. **结论边界：** 说明结果在哪些假设和数据范围内成立。

### 4.2 标题与摘要

- 标题直接包含模型、指标、算法或应用对象；
- 摘要通常按“问题 → 新方法/新结果 → 关键定量结果 → 化学意义”组织；
- 近期论文允许较长摘要，但不倾向把摘要写成完整 workflow 清单；
- “novel”“efficient”“promising”等词可以出现，但必须紧跟具体对象、比较基线或结果边界；
- 复杂方法论文会在摘要中点出主要数学结果，而不是只说“建立了一个框架”。

### 4.3 引言

- 先从图论、化学图、化学动力学或计算化学问题进入；
- 文献综述服务于当前定义或定理，不追求覆盖所有相邻领域；
- 通常在引言末尾给出清晰的 section roadmap 或贡献顺序；
- 化学联系需要落到分子、反应、拓扑指标或可计算对象上。

### 4.4 方法和结果

- 章节编号清楚，子章节标题承担导航作用；
- 定义、命题、定理、证明和公式集中出现；
- 图表用于验证具体结果，而不是同时承载背景、流程和结论；
- 计算实验会交代数据规模、参数、基线、硬件/数值方法和评价指标；
- 结论不会把数学界、数值验证界和化学应用界混成一个强结论。

### 4.5 版式与参考文献

- MATCH 自有页眉、卷期、DOI、作者和版权信息；
- 章节为编号形式，如 `1 Introduction`、`2 Main results`；
- 图题通常为 `Figure 1.`，表题为 `Table 1.`，公式持续编号；
- 参考文献使用编号列表，但条目按作者字母顺序排列；
- 文章可以较长，但页面密度来自公式、证明、图表和应用结果，而不是重复的解释段落。

## 5. 当前 `main.tex` 的事实审计

### 5.1 当前稿件的优点

- 已经明确把问题限定在固定拓扑、配位数相容构建单元、真旋转等价和刚体重构范围内；
- 已经区分 exhaustive、bounded prefix、sampled reconstruction audit 和 screening-stage DFT；
- Burnside、CASM、商图着色、稳定子、轨道和 Kabsch 都提供了可转化为 MATCH 形式化方法的材料；
- `ffc`、`fjh`、`vba` 形成了小空间、大空间和受限大空间的验证梯度；
- 27 个抽样指派的重构审计、CASM 交叉核对和固定组成 CO₂/N₂ 检查使化学计算联系不是空泛的。

### 5.2 当前稿件与 MATCH 的主要差距

| 优先级 | 当前事实 | 风险 | 修改方向 |
|---|---|---|---|
| 必须 | 使用 `sn-jnl` 和 `sn-mathphys-num` | 不符合 camera-ready 模板要求，可能在编辑初筛阶段被退回 | 下载/采用 MATCH 模板，迁移标题页、作者信息、章节、图表和参考文献环境 |
| 必须 | 使用 `\section*`、`\subsection*` | 文章导航不符合近期 MATCH 的编号章节习惯 | 改为编号章节，并采用 MATCH 风格的 `1 Introduction`、`2 Method/Results` 等结构 |
| 必须 | `main.bbl` 按正文首次引用顺序排列 | 不符合官方要求的字母排序 | 重建 BibTeX/bibliography 流程，按第一作者排序并逐项补齐书目信息 |
| 必须 | 稿件同时展开枚举、重构、DFT、GCMC 和多个应用接口 | 主贡献可能被读成泛化 MOF workflow | 将核心问题置于固定网络有序指派的非冗余枚举，其他内容作为验证层 |
| 强烈建议 | 轨道/稳定子/Burnside 逻辑主要位于叙述和算法中 | 数学正确性主张不够像 MATCH 的主结果 | 增加一个简短的 `Problem formulation and correctness` 小节，给出定义、命题和证明/验证边界 |
| 强烈建议 | 结果采用较多现代跨领域和筛选叙述 | 可能弱化与数学化学/计算化学的直接联系 | 用具体图规模、轨道计数、运行时间、覆盖范围和化学模型替代概括性宣传 |
| 强烈建议 | 章节标题偏 Nature/综合期刊风格，如 `Illustrative Screening-Stage DFT Example` | 与 MATCH 的直接标题习惯不完全一致 | 改成更短、更对象化的标题，如 `Single-point energy screening`、`CO2/N2 ordering application` |
| 可选 | 当前没有明确的复杂度/规模增长表述 | 算法贡献缺少可比较的计算解释 | 若代码和数据支持，补充搜索空间规模、轨道数、运行时间和瓶颈；不要无证据声称多项式复杂度 |
| 可选 | 仅给出代表性重构图和抽样审计 | 几何重构作为第二主贡献的定量支撑仍偏窄 | 可在补充材料或数据包中提供更完整的失败分类、输入清单和可复核日志 |

### 5.3 当前稿件中必须保留的科学边界

以下边界应在 MATCH 版本中继续保持：

- 只对固定输入网络和指定构建单元库进行非冗余枚举；
- proper-rotation equivalence 是取向兼容的默认规则，不声称覆盖所有 improper operation；
- `vba` 结果是按调度顺序得到的有界前缀，不是完整计数；
- 27/27 是抽样的 `ffc`/`fjh` 重构审计结果，不是所有设计空间的成功率；
- DFT 是刚性组装、未弛豫结构上的单点筛选分数；
- CO₂/N₂ 是固定协议下的筛选标签，不是最终吸附性能排名；
- 不把 sanitizable、topology-pass、clash-free 和 thermodynamic stability 混为同一层级；
- 不把不同枚举指派在后续弛豫后的几何独立性视为已被当前工作证明。

## 6. 推荐的 MATCH 目标结构

下面是实施阶段应采用的目标骨架。它不是新增科学结果，而是把现有内容改成 MATCH 更熟悉的证据顺序。

### 1 Introduction

用三段完成：

1. MOF 模块化和固定网络有序排布带来的组合问题；
2. 现有构建、识别、筛选方法与当前固定网络取商问题的边界；
3. 本文的算法贡献、重构贡献、验证层级和文章路线。

引言应避免把语言模型、生成模型、稳定性预测和吸附筛选都展开成平行背景；它们只需用于说明当前方法在计算流程中的位置。

### 2 Problem formulation and symmetry-constrained enumeration

建议包含：

- 输入：primitive quotient graph、degree-typed building-unit pools、composition counts；
- 配置空间 (X) 和 proper-rotation group \(\Gamma\)；
- group action、orbit、stabilizer、canonical representative；
- degree/composition pruning；
- subgroup-lattice scheduling 与 full-\(\Gamma\) canonical insertion；
- 一个短命题：在无运行上限且完整遍历时，canonical insertion 每个 proper-rotation orbit 保留一个代表元；
- Burnside 计数作为 exact closure 检查，而不是在 bounded mode 下的完整性保证。

### 3 Geometric reconstruction of periodic candidates

集中解释：

- assignment 如何提升回周期网络；
- node alignment 和 linker placement；
- Kabsch rotation、周期坐标处理、sanitization 和 topology audit；
- rigid-body/conformer scope；
- 输出是 pre-relaxation starting model。

### 4 Computational results

建议按证据强度排序：

1. `ffc` 的 exhaustive Burnside/CASM count checks；
2. `fjh` 的 larger-space matched bounded comparison；
3. `vba` 的 bounded-prefix behavior；
4. sampled `ffc`/`fjh` reconstruction audit；
5. runtime、space size 和 hardware/implementation details。

表格应优先呈现输入规模、naive assignments、orbit count、comparison count、time 和 run mode；不要让读者从正文长段落中拼接这些信息。

### 5 Chemical applications and downstream screening

将现有 DFT 和 CO₂/N₂ 内容放在一个明确的应用章节中：

- single-point energy screening；
- matched CO₂/N₂ ordering check；
- fixed-composition controls；
- screening-stage interpretation。

该章节必须在标题、首段、图注和结尾同步说明：这是对输出模型可进入下游计算的示范，不是稳定性或最终性能证明。

### 6 Discussion and limitations

讨论重点应是：

- 为什么取向兼容的非冗余枚举是 MOF 有序设计中的独立计算步骤；
- 与 substitutional enumeration、MOF construction/identification、screening workflows 的关系；
- 当前单一 primitive quotient graph、单一 embedding、rigid-body reconstruction、pre-relaxation CIF 和有限筛选样本的限制。

### 7 Conclusion

只保留三层结论：

1. 算法解决了什么固定网络枚举问题；
2. 现有证据在哪些范围内支持它；
3. 哪些问题必须通过 relaxation、larger adsorption studies 或更广泛数据验证。

在当前实施版本中，不再另设独立的 Methods 章：算法机制和实现细节保留在第 2 章，计算平台、基准体系和 CASM 对照放在第 3 章 `Results`，而 DFT 设置、CO₂/N₂ 模拟协议及重构审计放在第 4 章 `Reconstruction and Screening Validation`。随后安排 Discussion、Conclusion、Data availability、Acknowledgements、Author Contributions 和 Competing Interests，具体顺序以 MATCH 模板为准。

### 6.1 当前采用的章节结构

`paper/main_MATCH.tex` 已按以下结构完成第一轮迁移：

1. `Introduction`；
2. `Candidate Generation and Reconstruction`：无冗余枚举、几何重构、实现细节；
3. `Results`：基准与枚举协议、代数枚举核验；
4. `Reconstruction and Screening Validation`：重构审计、单点 DFT 筛选、CO₂/N₂ 排布核验；
5. `Discussion`：工作流意义与验证边界；
6. `Conclusion` 及投稿附加部分。

这一安排遵循“方法机制单独讲清楚，实验协议紧邻结果呈现”的原则，避免将 Python/硬件/软件版本、CASM、DFT 和 GCMC 设置混入算法定义，同时也避免用独立 Methods 章把实验结果与其评价条件割裂开。

## 7. 摘要与语言修改方向

### 7.1 摘要逻辑

推荐采用以下五句逻辑：

1. MOF 固定网络上的组分排布产生组合空间和对称冗余；
2. 现有流程没有同时给出取向兼容的非冗余指派和周期模型重构；
3. 本文提出 proper-rotation orbit enumeration + rigid-body Kabsch reconstruction；
4. 报告代表性网络中的 orbit-count/CASM/Burnside 结果和抽样重构审计；
5. 说明 DFT/CO₂/N₂ 仅作为下游筛选接口示例，并明确 pre-relaxation 边界。

摘要中应减少 `screening-ready`、`directly used`、`traceable route` 等可能被理解为成熟材料发现流程的表达，改用 `pre-relaxation starting models for downstream screening` 一类可核查表述。

### 7.2 术语表

| 推荐固定术语 | 不再交替使用的表达 |
|---|---|
| non-redundant enumeration | exhaustive generation / candidate generation（除非确实指对应模式） |
| proper-rotation group | orientation-preserving rotation group、rotation subgroup 混用 |
| quotient graph | reduced graph、primitive graph 随意替换 |
| orbit representative | retained assignment、canonical candidate 随意替换 |
| rigid-body reconstruction | geometric instantiation、assembly、placement 混用 |
| pre-relaxation starting model | screening-ready crystal、final candidate |
| single-point DFT screening score | DFT validation、energy ranking |
| fixed-composition comparison | pooled energy comparison |

## 8. 实施阶段验证清单

### 模板与编译

- [ ] 使用 MATCH 官方模板重新组织 `main.tex`；
- [ ] 标题、作者、单位、通讯作者和 received/date 字段符合模板；
- [ ] 所有章节和子章节编号；
- [ ] 图、表、算法和公式均有稳定编号；
- [ ] 无 undefined citation、reference、equation 或 rerun warning；
- [ ] 在 MATCH 模板下完成 PDF 编译。

### 参考文献

- [ ] 条目按第一作者字母顺序排列；
- [ ] 期刊名、卷、年份、页码/文章号和 DOI 与原始文献一致；
- [ ] MATCH 期刊名称统一写作 `MATCH Commun. Math. Comput. Chem.`；
- [ ] 正文引用顺序与编号列表重新生成并核对。

### 科学论证

- [ ] 明确固定网络问题的输入、输出和等价关系；
- [ ] 明确 exact/exhaustive 与 bounded/partial coverage；
- [ ] 明确 sampled reconstruction audit 的分母和范围；
- [ ] 明确单点 DFT、GCMC 标签和后续稳定性/性能验证之间的边界；
- [ ] 不新增实验，不把现有示例扩展成未验证的普适性结论。

### 视觉检查

- [ ] 首页和摘要版式；
- [ ] 算法跨页与行号；
- [ ] CASM 表格和 reconstruction audit 表格的可读性；
- [ ] 多面板重构图、DFT 图和 CO₂/N₂ 表格的图题/表题位置；
- [ ] 参考文献页的断行、缩进和 DOI；
- [ ] 最终 PDF 页面密度与 MATCH 近期论文相近。

## 9. 后续执行顺序

1. 取得并确认 MATCH 官方 LaTeX 模板，建立独立的投稿源文件，不覆盖现有 `sn-jnl` 版本。
2. 完成标题页、章节编号和参考文献格式迁移。
3. 先重排 Introduction、枚举方法和 reconstruction 方法，再调整 Results/Applications；不要先做逐句语言润色。
4. 加入简短的形式化问题定义和 correctness boundary，避免无证据增加复杂度定理。
5. 重新组织 DFT/CO₂/N₂ 为 chemical applications/downstream screening。
6. 编译、检查引用和交叉引用，完成 PDF 视觉检查。
7. 最后再做英文风格压缩，并同步更新 cover letter、Data availability、GitHub README 和 `codex/` 持久记忆。

## 10. 本轮不执行的事项

- 不修改 `paper/main.tex`；
- 不迁移或替换当前 LaTeX 模板；
- 不重新运行枚举、DFT 或 GCMC；
- 不新增定理、实验数据或材料性能结论；
- 不删除当前工作树中的既有修改。
