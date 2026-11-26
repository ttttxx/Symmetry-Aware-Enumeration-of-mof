# MATCH 全文术语校准：第一轮

日期：2026-11-10。审读对象：`paper/main_MATCH.tex`，含标题、摘要、正文、公式说明、算法、表格、图注和文末 Supplementary Information。

同日追查更新：见 `validation/MATCH_TERMINOLOGY_FOLLOWUP_2026-11-10.md`。PLD/LCD 已确认是解析器前两列映射颠倒，30/30 原始输出核验通过；下文 T01 的“尚未定位”是首次审查时状态，现已被追查结论取代。建议术语表包含多项原稿已正确使用的区分，不应整体理解为待修改错误清单。

作者决定及执行：topology recovery 保留原文；PLD/LCD 解析映射、当前论文数据与表格已修正；模板表已移至附录 A.1 并列全十种节点和固定连接体。其余建议尚不等同于作者认可的修改指令。

本轮逐段核对术语及其所指对象，提出后续逐句润色采用的规范；尚未改写主稿，也未完成语言终审。行号以本轮源码为准。源码 SHA-256：`05270b0be8a26f3200a339a1d7ed0a43117f6a2bca17e9bcc5ef9de31b9a5102`。

## 1. 先固定对象之间的关系

1. 固定输入：嵌入周期网络、所选重复晶胞、各边的连接体模板、按顶点度数分类的节点标签池。
2. 枚举对象：节点指派 `node assignment`，数学表示为顶点着色 `vertex coloring`，记为 $c\in X$。
3. 等价关系：所选保向嵌入网络操作诱导的顶点置换群 $\Gamma$ 对 $X$ 的作用。
4. 候选：一个 $\Gamma$-轨道；规范代表是该轨道中字典序最小的指派。
5. 构建输出：由选定规范代表和模板生成的周期原子起始模型；几何放置、连接修补和清理属于此步骤。
6. 后续评价：指定样本的连接/接触审计、GCMC 和单点 DFT。每种证据有各自的样本范围。

这六层在后续句子中保持明确。离散指派等价、原子模型等价、化学可实现性、弛豫后稳定性分别需要相应证据。

## 2. 建议采用的术语表

下表的“决定”是本轮建议规范，可据此讨论和敲定。保留有明确学科含义的不同词；同一对象不靠换术语制造句式变化。

| 概念 | 规范英文／中文 | 稿内变体或近邻词 | 决定与使用范围 |
|---|---|---|---|
| 构建单元总称 | building unit／构建单元 | building block, secondary building unit | 本文叙述统一 building unit；讨论原文术语时可保留 building block。SBU 仅用于确属次级构建单元的化学实体。 |
| 图上的点 | vertex／顶点 | node, site | 公式及图论操作用 vertex；说明周期位置用 site；node 留给节点构建单元或已定义的 node assignment。 |
| 节点构建单元 | node building unit／节点构建单元 | node unit, secondary building unit | 首次明确 node 的化学层级；后文 node template、node type 各指模板和类型。不能把所有有机节点都笼统叫金属 SBU。 |
| 节点类别 | node type／节点类型 | node label, node template | type 是类别，label 是类别的离散标识，template 是带坐标/构象的构建输入；三者关联但不互换。 |
| 全部顶点的类型指派 | node assignment／节点指派 | node-type assignment, labeling, ordering, pattern | 首次定义可用 node-type assignment；定义后以 node assignment 为主要名称。 |
| 指派的数学表示 | vertex coloring／顶点着色 | graph coloring, quotient-graph coloring | 在商图语境用 quotient-graph coloring。颜色表示节点类型标签，当前 $X$ 允许相邻顶点同色；无需 proper-coloring 约束。 |
| 具体空间排布 | ordered arrangement／有序排布 | ordering, ordered pattern | 引言和物理叙述可用；数学判等和计数优先 assignment。ordered 仅指规定周期内的确定性排布。 |
| 嵌入网络 | embedded periodic net／嵌入周期网络 | fixed MOF net, embedded net, topology | 先定义完整名称，再用 embedded net。包含坐标和晶格时不以 topology 代替。 |
| 网络的拓扑类别 | net topology／网络拓扑 | net, graph | ffc/fjh/vba 可以作为 net 名称；几何嵌入与拓扑类别分别描述。 |
| 周期商图 | quotient graph $Q$／商图 | quotient-graph space | $Q$ 是图，$X$ 是该图上的指派空间；写 assignment spaces on the quotient graphs。 |
| 周期选择 | chosen repeat cell／所选重复晶胞 | unit cell, supercell | repeat cell 为枚举用周期；仅确定与母胞关系时称 supercell。模拟晶胞另行说明。 |
| 数学相容约束 | degree-compatible／度数相容 | degree-typed, coordination-compatible | 约束为图度数与模板连接位点数相匹配。中文不用“配位数相容”代替“度数相容”。 |
| 节点/连接体的连接数 | connectivity; 3-connected node／连接度、三连接节点 | coordination number | 不将网络连接度自动解释为金属中心配位数。 |
| 标签集合 | degree-specific node-label pool／按度数分类的节点标签池 | admissible node labels, degree-typed pools | $\mathcal C_d$ 是标签集合；首次用 node-label pool，后文 label pool。admissible 只表示满足当前离散限制。 |
| 枚举空间 | assignment space $X$／指派空间 | coloring space, candidate space | $X$ 可称 coloring space；轨道去重前后区分 $X$ 和 $X/\Gamma$。 |
| 判等群 | orientation-preserving vertex-permutation group $\Gamma$／保向操作诱导的顶点置换群 | net symmetry group, rotation group | 首次完整限定来源；后文 permutation group 或 $\Gamma$。它是操作的置换像，可能含平移诱导的非平凡置换。 |
| 晶体学操作 | orientation-preserving operation／保向操作 | proper operation, proper rotation | 对含平移的仿射操作用 orientation-preserving；对 Kabsch 的矩阵用 proper rotation（真旋转）。 |
| 取向反转 | orientation-reversing operation／反向操作 | reflection | reflection 是其中一种；说明整个被排除操作类时用 orientation-reversing。 |
| 指派等价 | equivalent under $\Gamma$／在 $\Gamma$ 下等价 | symmetry-equivalent | 首次明确 acting group 后才可简写；不扩展到原子坐标或化学性质等价。 |
| 互异候选 | symmetry-distinct node assignments／对称互异节点指派 | distinct patterns, distinct orderings | 在数学方法中写 inequivalent under $\Gamma$ 最精确；结果的 orderings 须能追溯到轨道代表。 |
| 候选本体 | assignment orbit $\mathcal O_\Gamma(c)$／指派轨道 | candidate orbit, candidate | 一个轨道定义一个候选；不能把轨道数称为原始指派数。 |
| 另一类轨道 | vertex orbit $Hv$／顶点轨道 | H-orbit, orbit | $H$ 作用于 $V$ 的轨道，区别于 $\Gamma$ 作用于 $X$ 的指派轨道；上下文不明时写全称。 |
| 输出指派 | canonical representative／规范代表 | canonical ordering representative, canonical assignment, orbit representative | 涉及 $\kappa$ 用 canonical representative；一般 orbit representative 不自动包含字典序最小条件。 |
| 稳定子 | stabilizer subgroup／稳定子子群 | full stabilizer, stabilizer size | full 指相对于整个 $\Gamma$ 求得；基数写 stabilizer order。子群任务的阶数和完整稳定子阶数分别表述。 |
| $\rho$ | stabilizer fraction／稳定子比例 | symmetry index, relative size, exact $\rho$ | 定义为 $|\mathrm{Stab}_\Gamma(c)|/|\Gamma|$，不是晶体空间群阶数，也不能刻画全部对称结构。 |
| 固定点集合 | fixed-point set $\mathrm{Fix}(H)$／固定点集 | fixed-point domain, full-group/identity domain | 定义用 set；作为搜索域时可用 search domain。不引入数值固定点迭代的含义。 |
| 子群的偏序结构 | subgroup lattice／子群格 | lattice, descent | 与晶格 lattice 分开。descent 指当前方法命名，不意味着任务始终沿单一包含链运行。 |
| 精确分层 | exact-stabilizer-order stratum／精确稳定子阶数层 | exact-$\rho$ stratum, symmetry stratum, layer | 给定 $\Gamma$ 时 exact-$\rho$ 与 exact order 等价；正文用 stratum。low/medium/high 是区间分组，另称 $\rho$ bins。 |
| 输出限制 | per-stratum quota $B_s$／每层输出配额 | budget, per-$\rho$ quota | quota 是接受代表的数量上限；time limit、global cap 另列，不混写计算时间预算。 |
| 运行范围 | exhaustive enumeration; bounded enumeration／穷尽枚举、受限枚举 | complete enumeration, finite run | 所有 $X$ 本身有限；有限运行不能准确表达受配额/上限/超时限制。 |
| 算法操作 | canonicalization; canonical-prefix pruning／规范化、规范代表前缀剪枝 | canonical comparisons, prefix rejection | 区分完整指派规范化与部分指派剪枝。首次给 depth-first search (DFS)。 |
| 模型生成 | geometric construction／几何构建 | reconstruction, assembly | 全文主术语用 construction。assembly 可指具体片段组装操作。既有代码文件名、标签和仓库 URL 保持原标识。 |
| 节点定向 | Kabsch-based node alignment／基于 Kabsch 的节点配准 | Kabsch placement | Kabsch 求节点真旋转；后续平移完成放置。连接体轴向放置单独描述。 |
| 连接体 | linker; linker template／连接体、连接体模板 | edge template, edge building unit | 图论 edge 与化学 linker 分层；固定模板属于构建输入。 |
| 构建结果 | pre-relaxation periodic model／弛豫前周期模型 | crystal model, periodic structure, oriented periodic model | 首次说明 oriented 是模板取向处理；不能由此推断取向异构体枚举完备或热力学有序。 |
| 接触检查 | short-contact check／短接触检查 | clash check, nonbonded clashes | 证据表述用 short-contact check 并保留阈值；优化操作可写 clash reduction。 |
| 连接残差 | connection-site separation deviation／连接位点间距偏差 | closure error, bridge length | 首次定义 $|d_{XX}-1.485\,\mathrm{\AA}|$；bridge 的名义键表值、实际间距分别描述。 |
| 模型连接检查 | building-unit connectivity check／构建单元连接检查 | topology recovery | 按现有实现收窄至简化连接图比较，见第 3 节。 |
| 成分控制 | node composition／节点组成 | composition, node-label count vector | 指各节点类型计数；公式层可称 node-label count vector。固定连接体时可缩写 composition。 |
| DFT 比较 | single-point DFT energy comparison／单点 DFT 能量比较 | ordering energy, symmetry energy | 表述同组成、指定起始坐标的相对能量；不自动上升为平衡有序化能或稳定性排序。 |
| 孔描述符 | pore-limiting diameter (PLD); largest cavity diameter (LCD)／孔限制直径、最大空腔直径 | 稿内数值字段 | 名称保留常用含义；字段对应先核查，不能靠更名掩盖，见第 3 节。 |
| 吸附指标 | selectivity $S$; working capacity $WC$; working-capacity fraction $R$; adsorbent performance score $APS$ | regenerability 字段 | $R$ 按 $WC/q_{\mathrm{ads}}$ 定义，仅表示本压力摆动的释放比例。原始字段 regenerability 保留数据溯源关系。 |

## 3. 必须先处理的概念或证据问题

### T01：PLD/LCD 字段对应疑点（优先核查）

- 位置：595、615–635、849 行；数据包 `tables/ffc_co2n2_ordering_check.csv` 和 `targeted_zeopp_descriptors.csv`。
- 复核：30/30 个候选满足记录的 `pld_angstrom > lcd_angstrom`。例如首行分别为 27.67674 和 26.44087 Å，稿内组范围沿用了这些字段。
- Zeo++ 的孔几何区分 largest included sphere 与 largest free sphere；在相同模型及几何定义下，能沿通道移动的球也能被容纳，故自由球直径不应大于最大内含球直径。该不等式判断是根据几何定义作出的推论。[Zeo++ 作者源码仓库 README](https://github.com/mharanczyk/zeoplusplus/blob/master/README)
- 当前结论：存在字段映射、解析或定义/协议不一致的疑点。尚未定位原始 `.res` 和解析器，因此**尚不能断言两列必然互换**。
- 落地前提：核对原始 `-res` 输出的 Di、Df、Dif 与解析器列映射；确认后同步候选级数据、聚合范围、主表及文字中的范围。此项超出单纯文字校准，本轮未改数据。

### T02：topology recovery 大于当前检查的直接证据范围

- 位置：摘要 132 行，560–564、803 行，以及相应中文注释。
- 代码证据：`structure_validator.py:314–430` 根据 `bonds_all` 和构建单元归属组装连接图，把目标 `MultiGraph` 转成 `nx.Graph`，比较节点数、边数、度序列及图同构。此处未检查周期边的平移标记，且简化图丢失多重边信息。
- 建议规范：`the reconstructed building-unit connectivity graph matched the simplified target graph`；简写用 `building-unit connectivity check`。
- 可以保留本检查内的 graph recovery；全文“恢复目标周期拓扑”需要另有保留周期平移信息的检验支持。检查通过也不直接证明化学成键合理、弛豫稳定或合成可行。
- 该收窄需要同步摘要、审计结果、SI 协议和结论，不宜只改其中一句。

### T03：模板表的覆盖范围与筛选数据不一致

- 位置：表 2（479–499 行）标题称包含 benchmark、construction、screening 的模板，Results 开头又将所有算例指向该表。
- 数据证据：30 个吸附候选的 `composition` 共出现 10 种节点类型；表内只列出其中 5 种。筛选数据中额外出现 `3c_pyrrole_Ch`、`4c_In_1_Ch`、`4c_Ni_1_Ch`、`4c_V_1_Ch`、`4c_tetrathiafulvalene_Ch`。
- 建议：补齐模板和各自用途，或明确该表仅为某些枚举/构建/DFT 基准所用模板，另给吸附集清单。不要把不同样本池写成相同节点库。

### T04：保向约束与手性的论证需要限定

- 位置：229、419、457 行。
- 保向晶体学操作包括平移、旋转及螺旋操作，不能全部命名为纯 rotation；其线性部分行列式为 +1。[IUCr 晶体学操作说明](https://journals.iucr.org/j/issues/2018/05/00/in5013/index.html)
- `a reflection would reverse template handedness` 对手性模板才具有手性改变的化学意义；非手性模板不能由这句话推出出现不同对映体。
- 建议句：`We use orientation-preserving operations to define assignment equivalence, consistent with template placement by proper rotations.` 若补充手性说明，明确 `for a chiral template`。
- Kabsch 的取向保持与候选划分是方法选择及适用范围；由此无法自动证明所有得到的原子模型互不重合。当前候选定义保持在指派层。

### T05：相同 $\rho$ 的意义

- 位置：250–258、313–325、530–552、595、646、676 行。
- $\rho$ 只给稳定子阶数相对于 $|\Gamma|$ 的比例。相同阶数可对应不同稳定子子群；同一 $\rho$ 不意味着相同子群、相同空间群、相同原子对称性或相同孔结构。
- 建议保留匹配条件 `same node composition and stabilizer fraction`；646 行的 `symmetry index` 改为 `stabilizer fraction`。
- 图中的 low/medium/high 是 $\rho$ 区间分组，与 $\Omega_s$ 精确分层分别命名。

### T06：名义桥长、实际位点间距与两种 0.8 Å 判据

- 位置：448–452、563、806、835–840 行。
- 代码证据：`assembler.py:42–81` 的 closure 指标汇总 $|d-1.485|$ 的均方根及误差不超过 0.8 Å 的比例；融合操作使用实际位点间距严格小于 0.8 Å。
- 两个判据分别为 $d<0.8$ 和 $|d-1.485|\leq0.8$，数值相同不构成同一概念。
- 建议固定名称：`X-site fusion cutoff`、`connection-site separation deviation`、`nominal CIF bond-table length`。
- RMSD 在节点向量配准与间距偏差统计中均可使用，但每次写清参照对象。`bridge bond length` 不应让人误读为已优化或实测的化学键长。

## 4. 可直接执行的全文校准项

以下项目不改变算法或数值；术语规范敲定后统一落地并同步中文注释。

| 编号 | 位置（源码行） | 当前问题 | 建议处理 |
|---|---|---|---|
| T07 | 142、164、196、433 | building blocks/units，node unit/SBU 层级混杂 | 本文对象统一 building units；节点构建单元明确区分有机节点和金属 SBU。 |
| T08 | 195–227 | coloring 可能被理解为 proper coloring | 在首次数学定义处简短注明颜色是节点类型标签，$X$ 按度数池独立赋值。现有 2.1 增补曾回滚，不自动恢复旧句。 |
| T09 | 229、414–419 | permutation image 与 crystallographic operations 容易被当作同一群 | 保留 $\Phi$ 和 image 区分；$\Gamma$ 在本稿是有限置换群，不能直接叫完整晶体空间群。 |
| T10 | 274 | $\mathrm{Fix}(H)\cong\prod_P\mathcal C_{\deg(P)}$ 的结构类型未明 | 文字明确这是按每个顶点轨道的标签给出的自然一一对应；这些集合没有自动附带群结构。无需当作群同构解释。 |
| T11 | 261、285–289 | subgroup-lattice descent 可能暗示每步严格沿包含链下降 | 保留方法名；强调调度的是固定点任务，现有“相邻任务可不可比”说明应保留。 |
| T12 | 289、530、532 | DFS 首次展开缺少括号缩写；layer/stratum 混用 | 首次 depth-first search (DFS)；数学候选层统一 stratum，计算阶段用 task schedule。 |
| T13 | 312、315、324、332 | budget/quota 混用 | $B_s$ 统一 per-stratum quota，global output cap 与 task time limit 单独命名。 |
| T14 | 464–465 | quotient-graph spaces 暗示 $Q$ 自身是枚举空间 | 改为 assignment spaces on the quotient graphs。 |
| T15 | 465 中文 | degree-compatible 译成配位数相容 | 改回“度数相容”；金属配位数仅在明确金属中心/配位原子语境使用。 |
| T16 | 472、507、511、684 | assignment count 与 orbit count 不够稳定 | 分别写 raw assignment count $|X|$ 和 candidate/orbit count $|X/\Gamma|$；CASM 差异是作用群下的轨道数差异。 |
| T17 | 507–512、684 | 外部 CASM 比较可能被读作同判等条件验证 | 保留 $F$ 与 $\Gamma$ 的区别；ffc 数量一致只支持该计数一致，不能推出等价关系相同。 |
| T18 | 560、595、637、646、686 | ordering representative、canonical assignment、representative 混用 | 算法对象统一 canonical representative；物理排布可称 ordering，并明确模型源自该代表。 |
| T19 | 563–564 | closure audit、bridge lengths、fusion tolerance 过密交错 | 按 T06 分别定义名义值、间距残差和融合阈值；保留原统计数据。 |
| T20 | 643–646 | ordering energies 可能暗示有序化能或弛豫能 | 用 single-point energies of the constructed models 或 relative energies at fixed node composition。 |
| T21 | 676 | low/medium/high strata 与精确层混用 | 该图注使用 low/medium/high $\rho$ bins；明确端点约定。 |
| T22 | 803 | stabilizer size 与 order | 统一 stabilizer order；全稿涉及 subgroup 的 size 也按上下文检查。 |
| T23 | 812–840 | $s$ 同时表示稳定子阶数和晶胞扩张因子 | 建议扩张因子改为 $\lambda$，同步该协议的文字和阈值公式；枚举的 $s$ 保留。此为局部消歧，不改数值。 |
| T24 | 774 | proper pairwise joins 的修饰对象不清楚 | 写 pairwise joins that remain proper subgroups of $\Gamma$，明确 proper 修饰子群，不修饰 join 运算。 |
| T25 | 132、497、641、859–865 | DFT/PBE+D3/SCF 等首次定义层级 | 摘要与主文可分别独立定义 DFT；SI 首次展开 self-consistent field (SCF)。PBE+D3、基组、赝势保留准确软件名称及引用，不用泛称取代。 |
| T26 | 609、855 | $R$ 的正文名与 regenerability 字段名不同 | 数学定义清楚，可以保留 working-capacity fraction；报告数据映射，避免扩写成材料寿命、循环耐久性或普遍再生能力。 |
| T27 | 855 | product of means 与 mean of products | 保留主表 APS 为 mean S × mean WC；run-to-run SD 基于每次乘积。不能为统一名称而抹平两种统计量差异。 |
| T28 | 标题、132、164、182、586 | oriented construction/model 尚需操作性定义 | 建议保留标题，在方法首次说明涉及节点刚体定向和连接体轴向放置。不得暗示枚举了全部取向状态。 |

### 编辑残留

- 180 行：`\end{table}` 后出现未注释中文“我们将度数相容的节点指派表述为商图着色问题”，会进入正文；182 行对应中文注释含“在所选重复晶胞内，，并……”的缺文。应删除游离正文并恢复完整对照注释。
- 本轮未改该处，以保留本次术语审查所基于的作者工作稿；不应把此前“已编译通过”的记录当作当前工作稿没有编辑残留的证据。

## 5. 全文覆盖记录

| 范围 | 本轮检查重点 | 结果／对应编号 |
|---|---|---|
| 标题、摘要，105–134 行 | 枚举对象、保向、Kabsch 范围、审计证据、筛选样本 | T02、T04、T25、T28；候选完整覆盖未被错误加入摘要 |
| 引言及比较表，136–191 行 | building unit、占位类比、node/SBU、两项任务 | T07、T28；保留 occupation 与 MOF node assignment 的对象区分；标出编辑残留 |
| 2.1，194–258 行 | $Q/X/\Gamma/\Omega$、两类轨道、稳定子及 $\rho$ | T05、T08–T10；作用公式及轨道–稳定子关系相容 |
| 2.2 与算法 1，261–400 行 | fixed-point set、子群格、规范前缀、分层、配额、完备性条件 | T11–T13；$\mathrm{Fix}(H)$ 与精确稳定子层保持区分 |
| 2.3 与算法 2，402–460 行 | 仿射操作/置换像、真旋转、模板/连接/融合 | T04、T06、T09；不将连接体操作归入 Kabsch |
| 枚举结果及两张计数表，462–553 行 | 原始指派/轨道计数、群差异、受限覆盖 | T03、T12、T14–T17、T21；穷尽与受限输出边界保留 |
| 构建审计及模型图，556–589 行 | topology recovery、桥接残差、起始模型 | T02、T06、T18、T19、T28 |
| 吸附正文、公式、表，591–638 行 | 固定组成/比例、孔描述符、指标及均值定义 | T01、T03、T05、T18、T26、T27 |
| DFT 正文及图，640–678 行 | 单点能、组成、分层与排序边界 | T05、T20、T21、T25；21 ffc 与 32 fjh 不并入 72 例审计集合 |
| 结论及声明，680–706 行 | 结论术语与前文一致、弛豫/稳定性边界 | T02、T16、T18；数据公开状态仍需投稿前独立核验 |
| SI 对称群及生成元，727–769 行 | 占位结构空间群、仿射操作/置换、保向 | T04、T09；空间群符号不能直接充当原子模型对称群 |
| SI 配额证明，771–798 行 | 真子群/并生成、缓存、覆盖与前缀 | T13、T24；无超时和足够全局上限条件保留 |
| SI 审计与构建协议，800–844 行 | 简化连接图、两类阈值、扩张因子 | T02、T06、T22、T23 |
| SI 吸附协议，846–857 行 | 样本库、孔描述符、均值和 APS 口径 | T01、T03、T26、T27 |
| SI DFT 协议，859–865 行 | 计算术语、单点坐标条件 | T20、T25；closed-shell Cd/open-shell Cu 不能自动扩展为每个整个模型的电子态结论 |
| SI X 位点融合，867–896 行 | cutoff、periodic distance、端点合并、不同检查范围 | T06；0.8 Å 为操作值，未称为唯一最佳阈值 |

## 6. 下一阶段的逐句执行方式

1. 先落实 T01–T06 的证据边界和术语选择；如原始孔描述符尚未找回，显式保留待核验状态。
2. 将术语表落实到正文、表头、图注、算法说明及中文注释；保留公式标识符、文件名、原始字段名的溯源关系。
3. 再按标题→摘要→引言→方法→结果→结论→SI 推进逐句润色。给每句稳定编号，记录原句、修改句和理由；编号以实际句数为准，不机械凑成 200。
4. 每句检查术语所指、科学命题、证据范围、语法、指代与搭配；每段再检查中心命题、句间衔接和节奏。对需要数据或新论证的句子单独留问题，不用更顺的英文掩盖。
5. 英文确定后同步中文辅助注释；完成文本修改再编译并核对公式/引用、表宽和图注。

本轮检查：读取全文和相关实现；用 CSV 逐行计数验证 PLD/LCD 疑点及节点库差异；未修改主稿、数据、算法或图件，未运行 LaTeX 编译或科学计算。
