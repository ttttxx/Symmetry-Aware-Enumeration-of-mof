# 审稿意见修回审计与修改建议

日期：2026-08-24

最近同步：2026-08-31（已将固定组成 \texttt{ffc} CO2/N2 ordering 小实验补入 `paper/main.tex`；MTV-1–MTV-5 已用现有 CO2_N2 `batch_002_clean_pool_5000` 标签完成；主文改为 3 个 matched groups、共 12 个标签的紧凑摘要）

对象稿件：`paper/main.tex`

说明：本文件基于当前 `paper/main.tex` 与 `codex/` 项目记忆整理。它不是正式 rebuttal letter，而是修稿工作清单：哪些审稿意见已经部分处理，哪些仍需补强，以及建议如何修改正文、图件、数据可用性和回复口径。

## 一、当前总体判断

当前稿件已经完成以下重要修改：

1. 术语和缩写已明显补强：Abstract 中已展开 MOF、DFT、CASM；DFT 方法中已展开 HOMO-LUMO 和 SCF；正文新增了 `Terminology notes`，解释了 topological design inputs、symmetry-distinct assignments、proper-rotation group、Stab(c)、stratum、sanitizable periodic CIF 等。
2. 新颖性表述已降调：Introduction 已补充 ToBaCCo、MOFid/MOFkey、JACS 2025 free-energy ML screening、MOFGPT/MOFFlow，并把贡献限定为 fixed topological design space 内的 reconstruction-compatible symmetry reduction，而不是宣称生成全球全新的 MOF 或已验证稳定材料。
3. 2026-08-28 已完成 Fig. 2 / reconstruction audit 局部修订：Fig. 2 改为 `ffc/fjh/vba` 三个放大 panel；正文新增 reconstruction audit table，当前主文定量表只报告完整跑通的 sampled `ffc/fjh` 验证集合。
4. 2026-08-28 已完成标题降调、multivariate / digital reticular chemistry 补文献、Conclusion 新增和 Data availability 扩写：标题改为 `Screening-Ready MOF Candidate Models`；Introduction 新增 multivariate/digital reticular chemistry 段并引用 Lyu 2020、Zhang 2015、Ji 2020、Helal 2017、Canossa 2023；Discussion 后新增 `Conclusion`；Data availability 列出代表 CIF、`ffc/fjh` audit CIF、CSV 表格、manifest、validation summary、README、figure/table mapping 和 code-version notes。
5. 2026-08-28 已准备 Data availability package，且用户已确认对应 GitHub 数据包公开：`data_availability/mof-enumeration-reconstruction/` 本地副本包含 3 个 representative CIF、27 个 sampled `ffc/fjh` audit CIF、3 张 reconstruction CSV 表、`README.md`、`manifest.csv`、`validation_summary.md`、`figure_table_mapping.md` 和 `code_version.md`；正文 Data availability 已改为公开 GitHub URL。
6. 2026-08-28 第二轮正文修改已写入 `main.tex` 并编译通过（21 pages）：RCSR 输入边界与非 net classifier 说明（Methods）；related-work 对比表 `tab:related_work`（Introduction）；基准拓扑选取段（Results）；building-unit 清单 `tab:building_units` 与 Cd/Cu 说明（Methods/DFT）；colored-net + relaxation basin 第三限制段（Limitations）；Data availability 改为 definitive public URL；`reference.bib` 新增 `okeeffe2008rcsr`。
7. 2026-08-28 已完成全稿过强词与英文风格统一检查：`physically instantiable`、`prove/proof`、`fully correct`、`fully sanitized`、`promising/competitive/practical bridge` 等不再作为正文正向主张使用；保留的 `thermodynamic stability`、`synthesizability`、`relaxation` 仅用于限制、边界或未来验证语境；被改动英文段落均同步更新紧随其后的 `% 中文翻译：`。
8. 2026-08-31 已完成 R1.5/R3.5 的固定组成 ordering 小实验补文：`paper/main.tex` 新增 `Matched Ordering Check for CO2/N2 Screening` 小节和 `tab:ffc_co2n2_ordering`，使用 CO2_N2 `batch_002_clean_pool_5000` 中 3 个 matched comparison groups；每组固定同一 \texttt{ffc} composition、同一 `rho`、同一 symmetry stratum，并比较 4 个 symmetry-distinct ordering；Data availability package 新增 `tables/ffc_co2n2_ordering_check.csv`，共 12 个候选级标签。

如果目标是进一步提高修回成功率，仍有几类意见需要继续处理：

1. **数据包已按用户确认从待办中关闭**：当前工作区保留本地数据包 `data_availability/mof-enumeration-reconstruction/`，包含代表性 CIF、sampled `ffc/fjh` audit CIF、manifest、validation summary、summary CSV、README、figure/table mapping 和 code-version notes；正文现在使用公开 GitHub URL。
2. **Fig. 2 视觉问题已缓解，后续重点是文件可核查性**：主文图已放大并拆成三个 panel；现在需要保证 Fig. 2 对应 CIF 和 Table 相关 CSV 能被审稿人直接下载核对。
3. **DFT 标题风险已处理，正文仍需保持边界一致**：标题已去掉 “Physically Instantiable”；正文继续把 single-point DFT 限定为 screening-stage score，不写稳定性、可合成性或 relaxed identity。
4. **R2 的三个技术问题已补入正文并编译通过**：proper-rotation group 的物理合理性、Fix(H) overlap / double-counting、interpenetrating / multiple-cell / non-simply connected graphs 已分别加入 `Organizing Non-Redundant Candidate Generation`、Algorithm 1 后说明段和 `Current Limits of the Validation`。
5. **Conclusion 缺失已处理**：正文已在 Discussion/Limitations 后、Methods 前新增 `Conclusion`，三段分别对应方法贡献、证据边界和后续 relaxation/property screening。
6. **英文风格与过强词已完成一轮全稿收紧**：Abstract、Introduction、Results、Discussion、Conclusion、Methods 和 Data availability 均已按中性边界重写。
7. **晶体工程实证已补入正文**：方案 C 已收窄为 existing-label 小实验，不新增大规模吸附研究；主文报告 3 个同组成、同 \(\rho\)、同 symmetry stratum 的 matched groups，显示不同 ordering 的孔结构描述符和 CO2/N2 screening score 存在差异。

## 二、逐条意见跟踪表

列说明：

- **为什么这么修改**：回应审稿人关切背后的逻辑，不是操作步骤。
- **已做正文修改**：写明 `main.tex` 位置 + 中文一句简述；代码/数据锚点见括号。
- **下一步 TODO**：正文已完成写 **无**；实验/push/可选 polish 写具体 ID（见第八节）或一句话待办。

| ID         | 审稿人意见核心                                                                             | 为什么这么修改                                                                                                         | 已做修改                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | 下一步 TODO                                                                                                                                 | 优先级 |
| ---------- | ------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- | ------ |
| R1.1       | Fig. 2 显示配体叠加；“sanitizable periodic models” 表述不真实；缺少 CIF 证明             | 审稿人需要可下载、可核查的 CIF 与定量验证，而非仅靠放大图；同时不能把 pre-relaxation starting model 说成已验证稳定材料 | **Fig. 2**（`fig:reconstruction_gallery`）：三 panel 放大 ffc/fjh/vba；图注写明 visual examples 不作稳定性证据（中图仅作可视化，vba 不进定量表）<br>**Table** `tab:reconstruction_audit`：ffc 12/12、fjh 15/15，合计 27/27 sanitizable/topology-pass/clash-free（定量审计仅限完整跑通的 ffc/fjh 抽样集）<br>**Data availability**：已改为公开 GitHub URL，并列出 package 内容；代码锚点：本地包 `data_availability/mof-enumeration-reconstruction/`；2026-08-28 审计输出 `results/reconstruction_summary_split/partial_2026-08-28/` | 无 | 高     |
| R1.2       | rotations 不能普遍解决 colored net 区分；适用范围不清；需分析 RCSR/EPINET/ToposPro         | 审稿人混淆了「跨库 net 分类」与「固定 net 内 assignment 去冗余」；需明确 RCSR 是输入来源，本文不是 net classifier      | **Methods §Computational Details**：`workflow starts from…RCSR…not a net classifier`（基准 net 来自 RCSR embedding；方法只在固定输入内做 building-block 指派去冗余，不判断两个 net 是否同色同构）<br>**Limitations 第三段**：colored-net 分类须在 workflow 前完成旋转不能替代跨库/跨 embedding 的 net 等价判别）<br>代码锚点：`data/mof/{topology}.cif` → graph（`experiment/casm_compare_README.md`）；引用 `okeeffe2008rcsr`                                                                                                      | 无                                                                                                                                          | 高     |
| R1.3       | single-point DFT 不能得到 physically instantiable candidates；需要晶胞和坐标优化           | 原标题/DFT 表述易被理解为已物理可合成或热力学稳定；需把证据边界收窄到 screening-stage                                  | **标题**：`Screening-Ready MOF Candidate Models`（候选模型可进入筛选，非已物理实例化）<br>**DFT Methods + Results + Conclusion**：single-point screening score；final ranking requires relaxation（DFT 是组装几何上的筛选分数，不是稳定性/可合成性证明）                                                                                                                                                                                                                                                                                      | 无                                                                                                                                          | 无     |
| R1.4       | 概念和语言不严格，术语难懂                                                                 | 化学/晶体学读者需要可读的符号定义，避免 announcement-style 长句                                                        | **Abstract / Terminology notes / 各节 `%中文翻译`**：MOF、DFT、CASM、HOMO-LUMO、SCF、Stab(c)、Γ、stratum、sanitizable 等已展开或定义；全稿已完成一轮过强词和英文风格统一检查，长句和防守式表达已压缩                                                                                                                                                                                                                                                                                                                        | 无                                                                                      | 中     |
| R1.5       | 晶体工程影响不清；为何选 ffc/fjh/vba；未连接 multivariate/digital reticular chemistry      | 需解释 benchmark 拓扑不是随机选取，且与 multivariate ordered assignment 的应用场景相连；需用同组成不同 ordering 的性质差异作最小实证           | **Introduction**：multivariate/digital reticular chemistry 段（功能取决于 building block 在晶体学位置上的排布，不仅是组分清单）<br>**Results 开头**：ffc=7 节点穷举核对，fjh=28 节点双固定组成 DFT，vba=bounded 压力测试（三拓扑分工明确）<br>**Results §Matched Ordering Check**：新增 3 个 \texttt{ffc} matched groups，共 12 个 CO2/N2 screening labels；每组固定 composition、\(\rho\)、symmetry stratum，只比较 symmetry-distinct ordering；APS 组内范围为 0.588–0.760、0.595–0.817、0.626–0.839<br>数据锚点：`data_availability/mof-enumeration-reconstruction/tables/ffc_co2n2_ordering_check.csv` | 无 | 中     |
| R1.6       | 英文质量低，图太小                                                                         | 图件可读性与 prose 质量直接影响 R1 整体印象                                                                            | **Fig. 2** 已改为三放大 panel 并编译检查（主文图件尺寸问题已缓解）<br>**Data availability** 已扩写 package 清单并改为公开 GitHub URL<br>**Language pass**：Abstract、Introduction、Results、Discussion、Conclusion、Methods 已按中性、精炼、非宣传式语气收紧，并同步中文翻译                                                                                                                                                                                                                                                                       | 无                                                                                 | 中     |
| R2.1       | Abstract 首次出现 MOF、DFT 需全称                                                          | 期刊格式要求                                                                                                           | **Abstract** 首次出现处已写全称                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | 无                                                                                                                                          | 无     |
| R2.2       | Results 句子中 DFT 首次出现需全称                                                          | 同上                                                                                                                   | **Results §Illustrative Screening-Stage DFT Example** 首次出现已写全称                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | 无                                                                                                                                          | 无     |
| R2.3       | DFT section 中 HOMO-LUMO、SCP/SCF 首次出现需全称                                           | 同上；reviewer 的 SCP 应为 SCF                                                                                         | **DFT Methods**：HOMO-LUMO、SCF 已展开（方法段缩写已规范）                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            | 无                                                                                                                                          | 无     |
| R2.4       | 为什么只用 proper-rotation group 而不用 full crystallographic point group？                | 这是 oriented reconstruction 的保守默认，不是 universal physical equivalence                                           | **Organizing Non-Redundant Candidate Generation**（Γ 定义后）：proper rotations 为保守重构选择；mirror/inversion 对取向多原子 BU 不一定等价（默认只用真旋转取商，避免重构后不等价的指派被提前合并）                                                                                                                                                                                                                                                                                                                                                  | 无                                                                                                                                          | 无     |
| R2.5       | 不同子群 fixed-point spaces 是否重叠？算法如何避免 double-counting？                       | 审稿人担心子群任务重叠导致重复计数                                                                                     | **Algorithm 1 后说明段**：Fix(H) 可重叠；subgroup lattice 仅调度；最终由 full-Γ canonicalization + Ω/Trie 去重（重叠子群任务不会产出重复候选）                                                                                                                                                                                                                                                                                                                                                                                                      | 无                                                                                                                                          | 无     |
| R2.6       | 是否适用于 interpenetrating structures、multiple unit cells、non-simply connected graphs？ | 需明确当前实现边界，避免过度承诺                                                                                       | **Current Limits of the Validation**：默认 single primitive quotient graph + one embedded net；互穿/超胞需显式 multi-net 或 supercell input（当前不自动生成互穿或多 net；非单连通主要影响重构验证）                                                                                                                                                                                                                                                                                                                                                   | 无                                                                                                                                          | 无     |
| R3.1       | 与 Rosen/Snurr/Kulik 等已有工作的区别不清                                                  | 需把本文定位在 candidate preparation，而非替代 ToBaCCo/MOFid/ML screening                                              | **Introduction**：`We do not replace ToBaCCo, MOFid/MOFkey…` + **Table** `tab:related_work` 四行对比 assembly / identifier / ML / 本文（本文在下游流程前做固定 net 上的非冗余有序指派与 oriented reconstruction）                                                                                                                                                                                                                                                                                                                          | 无                                                                                                                                          | 无     |
| R3.2       | 无输出 CIF，无法验证                                                                       | 审稿人无法独立核查结构与表格                                                                                           | **Data availability**：已给出公开 GitHub URL，并列 representative CIF、27 个 ffc/fjh audit CIF、CSV、README、manifest（正文已描述完整数据包结构）<br>本地包副本：`data_availability/mof-enumeration-reconstruction/`                                                                                                                                                                                                                                                                                                                                  | 无                                                                                                                            | 高     |
| R3.3       | single-point DFT 会保留不合理结构；relaxation 后候选可能合并                               | 需承认未 relaxation 的 interpretive 风险；若有实验只能写收敛且 topology-retained 的子集                                | **正文部分 — Limitations 第三段**：distinct assignments 可能弛豫到同一 basin；单点 DFT 不能证明每个 symmetry class 对应互异 relaxed material（弛豫后名义不同候选可能合并）<br>**实验部分**：正文未写 “we performed relaxation”                                                                                                                                                                                                                                                                                                               | **EXP-1/2**：`dft_validation/run_relaxation_validation.py`；当前 `ffc_comp_a_11` 为 `harvested_partial_restart`，不能作正面证据 | 高     |
| R3.4       | 只生成候选，没有稳定性/可合成性；不知道输出可靠性                                          | 需提供 output reliability 的可核查指标，并明确不证明稳定性/可合成性                                                    | **Table** `tab:reconstruction_audit` + 表注：attempted/assembled/sanitizable/topology-pass/clash-free（指标衡量 starting-model generation 可靠性，不是稳定性）<br>Data availability 已给出公开 GitHub URL；本地 README/manifest/validation summary 已备                                                                                                                                                                                                                                                                                         | **AUDIT-1**（可选，仅当 vba 需进同表；当前 vba 10 个 coloring 仅 1/10 完整四项验证）                                | 高     |
| R3.5       | promised relevance to crystal engineering does not appear                                  | 需把方法接到 multivariate ordered assignment 与后续 screening 接口，并用同组成 ordering 对比显示 ordering 可影响筛选性质                               | **Introduction** multivariate 段 + **Conclusion** 后续 relaxation/GCMC/ML 用途（方法是 fixed-net 上有序指派的候选准备，不是稳定材料发现）<br>**Results §Matched Ordering Check**：3 个同组成 matched groups 显示 PLD/LCD、CO2 working capacity、selectivity、regenerability 和 APS 在 ordering 间有差异；口径限定为 screening-stage relevance，不写最终吸附性能排序                                                                                                                                                                                                                                              | 无                                                                 | 中     |
| R3 minor 1 | CASM 是什么                                                                                | 缩写需解释                                                                                                             | **Abstract + Results**：CASM 首次出现已展开（Clusters Approach to Statistical Mechanics 已给出）                                                                                                                                                                                                                                                                                                                                                                                                                                                      | 无                                                                                                                                          | 无     |
| R3 minor 2 | Fig. 2 unit cells 太小                                                                     | 主文图需足够大以便审查结构                                                                                             | **Fig. 2** 三放大 panel                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | **SI-1**（可选）：SI 全尺寸结构图 + 对应 CIF                                                                                          | 低     |
| R3 minor 3 | 本文考虑了哪些过渡金属                                                                     | 不能只写 “transition metals”                                                                                         | **Methods §Building-Unit Inventory** + **Table** `tab:building_units`：主文 DFT 以 **Cd** 为主（闭壳，`USE_CLOSED_SHELL_VARIANT=True`），**Cu** 为 open-shell probe 与 vba 示例池（明确列出 Cd/Cu 及实际 nodes/edges）<br>代码锚点：`dft_validation/config.py`；`experiment/geometric_reconstruction_summary.py` CONFIGS                                                                                                                                                                                                   | 无                                                                                                                                          | 无     |
| R3 minor 4 | 无 Conclusion                                                                              | 稿件结构不完整                                                                                                         | **Conclusion** 三段：方法贡献 / 证据边界 / 后续 relaxation·screening（结论已补，且不夸大为 discover stable MOFs）                                                                                                                                                                                                                                                                                                                                                                                                                                    | 无                                                                                                                                          | 无     |

## 三、建议直接修改 `main.tex` 的位置

### 1. 标题

原标题：

`From Topological Design to Physically Instantiable MOF Candidates: Symmetry-Aware Enumeration and Reconstruction`

已改为更保守、更抗审稿质疑的标题：

`From Topological Design to Screening-Ready MOF Candidate Models: Symmetry-Aware Enumeration and Reconstruction`

原因：`Physically Instantiable` 被审稿人抓住后，会被理解为“已经物理稳定 / 可合成 / relaxation 后仍成立”。当前证据更适合支持 “screening-ready candidate models” 或 “reconstructable candidate models”。

### 2. Abstract

当前 Abstract 已经降调并与 reconstruction audit 口径同步：主文不再写 “sanitizable periodic CIFs across five topologies”，而是写 completed sampled `ffc/fjh` audit gives `27/27` sanitizable, topology-consistent, and clash-free periodic starting models。仍需保持以下解释边界：

- `sanitizable` 只表示格式和重载通过；
- `topology-pass`、`clash-free` 是更强的验证；
- DFT 只说明进入 screening-stage evaluation，不说明稳定性。

如后续需要继续压缩 Abstract 或 rebuttal，可使用这一句作为边界句：

`These CIFs should be interpreted as validated starting models with recorded reconstruction diagnostics, rather than as relaxed or synthesis-confirmed MOFs.`

### 3. Fig. 2 / reconstruction validation

当前 Fig. 2 / reconstruction validation 已按“两层证据”处理：

1. 主文保留 3 个清晰放大的代表结构：`ffc`、`fjh`、`vba`；
2. 定量表只列完整跑通的 sampled `ffc/fjh` audit metrics：
   - attempted
   - assembly success
   - sanitizable/reloadable CIF
   - topology-pass
   - clash-free before relaxation
   - representative CIF path

当前 2026-08-28 `ffc/fjh` 重跑结果可作为主文表格依据：

- `ffc`: attempted `12`, assembly `12/12`, sanitizable `12/12`, topology-pass `12/12`, clash-free `12/12`
- `fjh`: attempted `15`, assembly `15/15`, sanitizable `15/15`, topology-pass `15/15`, clash-free `15/15`
- combined current sampled audit: attempted `27`, assembly `27/27`, sanitizable `27/27`, topology-pass `27/27`, clash-free `27/27`

注意：7 月旧结果中较低的 topology/clash 通过率来自 assembler/validator 修复前的审计，不应再作为当前主文表格数字。2026-08-28 另做过 `vba10` 诊断尝试，10 个 `vba` coloring 中只有 1 个完整走完 assembly/sanitize/topology/clash 四项验证，另外 9 个在 native 层退出，未形成完整 metadata row。因此主文不要把 `vba` 写成已完成同等定量审计；当前稳妥口径是“`vba` 为较大 bounded-regime visual example，定量 audit 限于完整完成的 `ffc/fjh` sampled validation set”。

### 4. Proper-rotation group 解释

已在 `Organizing Non-Redundant Candidate Generation` 中、`Γ` 定义段之后加入该说明段：

`This proper-rotation restriction is a conservative reconstruction choice, not a claim that the full crystallographic point group is never appropriate. Improper operations can be valid equivalences for some building-unit libraries. In the present setting, however, building units are polyatomic objects with fixed connection vectors and input conformers. A mirror or inversion can change their handedness or local binding orientation after placement. For this reason, the default workflow quotients by orientation-preserving rotations before reconstruction. If a library is known to be invariant under a specific improper operation, that operation can be added as a system-specific equivalence rule.`

当前正文口径：这不是说 full point group 错，而是本文默认避免把重构后不一定等价的结构提前合并。该段不需要新增参考文献，因为它解释的是本文算法默认等价关系。

### 5. Fix(H) overlap / double-counting

已在 Algorithm 1 后、`Terminology notes` 前加入该说明段：

`Fixed-point spaces from different subgroups can overlap. A coloring with a larger stabilizer can belong to Fix(H1) and Fix(H2) at the same time. The subgroup lattice is used for scheduling, not as a disjoint decomposition of the design space. When a complete coloring is produced, it is canonicalized over the full proper-rotation group Γ. Only the canonical representative is inserted into the global set Ω and the Trie. Overlapping subgroup tasks may revisit the same candidate, but they do not create duplicate outputs. The orbit-partition fingerprint removes only tasks with identical vertex-orbit constraints. Final uniqueness comes from canonical insertion.`

这个段落可直接用于回答 R2.5。它明确区分了“子群任务调度”和“最终输出去重”。

### 6. Interpenetration / multiple cells / non-simply connected graphs

已在 `Current Limits of the Validation` 中补入该限制段：

`The current implementation also assumes a single primitive quotient graph and one embedded net. It does not automatically generate interpenetrated frameworks, multiple independent nets in the same cell, or supercell colorings. These cases require the corresponding multi-net or supercell quotient graph as explicit input, including the relative registry between net copies. Non-simply connected cycles do not change the enumeration principle, because the group action is still defined on quotient-graph vertices. The practical limitation is in reconstruction, where periodic consistency, interpenetration contacts, and long-range geometric checks need extra validation.`

这个段落能回答 R2.6，同时不做过度承诺。当前没有新增 interpenetrated MOF 参考文献，因为本轮只是说明当前输入模型和实现边界，不扩展到互穿结构综述。

### 7. Multivariate / digital reticular chemistry

#### 7.0 这两个概念是什么

**Multivariate reticular chemistry（多变量网状化学）**

- 研究对象：在**同一拓扑网络**上，同时使用**多种** metal node / organic linker，且它们的**空间排布（ordering）** 是有意设计的 MOF，而不是“有哪些组分”的简单混合。
- 核心问题：功能（催化、吸附选择性、光学响应等）往往不只取决于化学式，还取决于不同 building block **落在哪些对称不等价位点**、以及局部 **序列/排列** 如何分布。
- 典型文献：Zhang et al., JACS 2015（multivariate MOF-177）；Helal et al., Natl. Sci. Rev. 2017；Ji et al., Science 2020（metal sequencing）；Canossa et al., Nat. Rev. Mater. 2023（multivariate reticular sequences）。

**Digital reticular chemistry（数字网状化学）**

- 由 Lyu et al., Chem 2020 提出的更大框架：把网状化学里的 **设计空间搜索、自动化生成、数据库化与高通量筛选** 做成可计算、可复现的工作流。
- 它不是一个具体 MOF 案例，而是说：reticular design  increasingly 是 **大规模 combinatorial search + automation** 问题；需要工具把“可能的组分排布”变成 **可检查、可下游评估** 的结构模型。

#### 7.1 本文与它的关系（不是纯口径收紧）

| 层面 | 是否已有联系 | 说明 |
| --- | --- | --- |
| **算法/形式化** | ✅ 有 | 商图 **coloring** = 给对称不等价位点分配不同 building unit；即 fixed-net 上的 **ordered multivariate assignment** 枚举 + 对称去冗余 + 取向重构 |
| **现有 benchmark** | ✅ 有 | `ffc` 使用多类 node/linker；主文 DFT 是 fixed-composition symmetry screening，新增 CO2/N2 小实验则按 fixed-composition ordering check 组织 |
| **crystal-engineering 实证** | ✅ 已补 | 现有 CO2_N2 `batch_002_clean_pool_5000` 标签显示：3 个 matched \texttt{ffc} comparison groups（共 12 个候选）在固定 composition、\(\rho\)、symmetry stratum 后，不同 ordering 仍具有不同孔结构描述符和 CO2/N2 screening score |

正文已写入口径（Introduction 82–83 行、Conclusion 431–432 行）：

- multivariate / digital reticular chemistry 说明“功能不仅取决于有哪些 building blocks，还取决于它们如何排布在 crystallographic positions 和 local sequences 上”；
- 本文 workflow 提供 fixed net 上的 non-redundant, inspectable starting models。

**边界（修回时必须保持）**

- 可说：coloring ≡ ordered assignment；workflow 为 multivariate ordered assignments 的 **candidate preparation**；
- 不可说：已解决 experimental sequencing、已证明 sequence-dependent synthesizability、或 digital reticular 全空间自动化平台。

Rebuttal 可用压缩句（MTV 实验完成前；已不再作为当前版本主口径）：

`For multivariate reticular chemistry, the practical question often includes both building-block identity and ordering over symmetry-related positions. Our quotient-graph coloring formalism enumerates such non-redundant ordered assignments on a fixed net and reconstructs them into explicit starting models, which can then be passed to relaxation, adsorption simulation, or active-learning screening.`

Rebuttal 可用压缩句（**MTV 实验完成后，当前版本采用**）：

`On the ffc net, we further compared three matched fixed-composition groups using existing CO2/N2 GCMC-derived screening labels. Each group contained four symmetry-distinct orderings with the same composition, rho value, and symmetry stratum. The ordered assignments differed in pore descriptors and derived screening scores. We present this as an ordering-sensitivity check for candidate preparation, not as a final adsorption-performance ranking.`

#### 7.2 方案 C（已完成：MTV-1–MTV-5）

完成方式：没有新增大规模吸附研究，而是复用 `/Users/jing/PyCharmMiscProject/CO2_N2` 中已完成的 `batch_002_clean_pool_5000` 标签。主文选用 3 个 matched `ffc` comparison groups；每组固定 composition、\(\rho\)、symmetry stratum，并包含 4 个 symmetry-distinct ordering，共 12 个 GCMC-derived labels。主文报告 PLD/LCD、CO2 working capacity、CO2/N2 selectivity 和 APS 的组内范围；CSV 保留 regenerability 和完整候选级明细。

主文位置：`Results` 中新增 `Matched Ordering Check for CO2/N2 Screening` 小节和 `tab:ffc_co2n2_ordering`。Methods 新增 Zeo++/RASPA2/CO2-N2 screening label 定义和参考文献，Data availability package 新增 `tables/ffc_co2n2_ordering_check.csv`。

不画图理由：该检查用于回应 crystal-engineering relevance 的接口意义，不是论文主结果。主文保留 3 行范围表，完整 12 行候选级数据放入 CSV，避免把读者注意力从 symmetry-aware enumeration / reconstruction 主线转移到吸附筛选项目。

边界：该结果只说明 fixed-composition ordering 可改变筛选阶段描述符和派生分数；不作为 final adsorption-performance ranking、thermodynamic stability 或 synthesizability 证据。

#### 7.3 正文已引文献

- Lyu et al., Chem 2020, digital reticular chemistry
- Zhang et al., JACS 2015, multivariate MOF-177
- Ji et al., Science 2020, sequencing of metals
- Helal et al., National Science Review 2017, multivariate MOFs
- Canossa et al., Nature Reviews Materials 2023, sequences in multivariate reticular structures

### 8. DFT 与 relaxation

建议在 DFT 结果段和 limitations 中继续收紧：

- 避免 “physically instantiable”；
- 明确 single-point DFT 不证明稳定性；
- 明确 relaxation 后不同 candidates 可能并入同一 basin；
- 如果不做新计算，说明本文只给出 screening-stage retention test；
- 如果能补实验，优先选 Fig. 2 中被质疑的结构和 DFT minima 做少量 cell+atomic relaxation，作为补充表或 SI。

建议回复口径：

`We agree that single-point DFT cannot establish thermodynamic stability or synthesizability. We have revised the manuscript to describe these calculations only as screening-stage scores on as-assembled geometries and to state that full cell and coordinate relaxation is required before thermodynamic ranking.`

如果没有 relaxation 结果，不能写 “we performed relaxation”。

### 9. Conclusion

已新增 `\section*{Conclusion}`，放在 Discussion/Limitations 后、Methods 前。内容控制在 3 段：

1. 总结方法贡献：proper-rotation orbit enumeration + oriented reconstruction；
2. 总结证据边界：Burnside/CASM 检验、代表性 CIF、reconstruction diagnostics、DFT only screening-stage；
3. 说明用途和下一步：multivariate ordered assignments、relaxation、GCMC/DFT/ML active learning。

不要在 Conclusion 里重新夸大为 “discover stable MOFs”。

## 四、需要在单独 Git 数据仓库补充的数据文件与说明

为了回应 “CIF files were not provided” 这一数据可核查性问题，需要准备以下材料，并与正文 Data availability 保持一致：

本地已准备的数据包位于：

`data_availability/mof-enumeration-reconstruction/`

公开仓库中建议保持或等价保留以下结构：

1. `cifs/representative/`：主文 Fig. 2 每个 panel 对应的 CIF。
2. `cifs/reconstruction_audit_ffc_fjh/`：2026-08-28 sampled `ffc/fjh` audit 的 27 个 CIF。
3. `tables/reconstruction_records.csv`：每个 CIF 的 topology、candidate id、rho、stabilizer order、reloadable/sanitizable、topology_pass、clash_pass、warning message。
4. `tables/reconstruction_by_topology.csv` 和 `tables/reconstruction_overall.csv`：按 topology 和 overall 汇总 attempted / sanitizable / topology-pass / clash-free。
5. `README.md`：说明 sanitizable、topology-pass、clash-free 的定义，以及这些结构是 pre-relaxation starting models。
6. `figure_table_mapping.md` 或同等文件：说明正文 Fig. 2 和 reconstruction audit table 对应哪些 CIF/CSV。
7. `manifest.csv`：列出数据包文件、来源路径和用途。
8. `validation_summary.md`：用人可读方式总结 sampled `ffc/fjh` audit 的 27/27 结果和证据边界。
9. `code_version.md`：记录代码来源、commit / release tag，以及数据包准备时的本地修改状态。
10. 如果未来有收敛且 topology-retained 的 relaxation，再补充 `relaxed_subset/`；当前不要伪装成已有 relaxation 结果。

## 五、建议 rebuttal 总策略

对 R1/R3 的语气建议采取“承认边界 + 具体修稿 + 不硬辩”的方式：

- 承认原稿把 single-point DFT 和 “physically instantiable” 说得过强；
- 承认 Fig. 2 和数据可用性不足；
- 强调已将论文目标收窄为 candidate generation / reconstruction / screening preparation；
- 用 CIF、manifest、reconstruction audit、图件放大来提供可核查证据；
- 对不可完成的大规模 relaxation，不要以工作量为理由拒绝，而是说本文结论不依赖 relaxed stability ranking，并把 relaxation 放入后续验证协议。

对 R2 的语气可以更直接：

- 感谢 technical questions；
- 每个问题都在 Methods/Discussion 中新增一段；
- 明确 proper rotations 不是 universal equivalence；
- 明确 Fix(H) overlap 由 canonicalization + global Trie 去重；
- 明确 interpenetration/supercell 不是当前默认输入，但可通过显式 multi-net/supercell quotient graph 扩展。

## 六、执行顺序（2026-08-28 更新）

### 已完成（正文此前已有）

1. 改标题，去掉 `Physically Instantiable`
2. proper-rotation、Fix(H) overlap、interpenetration/multiple-cell 三段技术说明
3. Conclusion、multivariate/digital reticular chemistry
4. Fig. 2 三 panel + ffc/fjh reconstruction audit 表

### 本轮：正文修改已写入 `main.tex`（2026-08-28）

5. **§7.1** RCSR 输入边界 + 非 net classifier（Methods）— 已完成
6. **§7.2** related-work 对比表 `tab:related_work`（Introduction）— 已完成
7. **§7.3** ffc/fjh/vba 选取理由（Results 开头）— 已完成（hwx/nbo 已从正文移除）
8. **§7.4** building-unit / Cd·Cu 清单 `tab:building_units`（Methods）— 已完成
9. **§7.5** colored-net + relaxation basin 边界（Limitations）— 已完成
10. **§7.6** Data availability definitive public URL 措辞 — 已完成
11. **`reference.bib`** 新增 `okeeffe2008rcsr` — 已完成
12. **`main.tex` 编译通过**（21 pages）— 已完成

### 本轮已关闭 / 仍保留

- GitHub 数据包 push 与链接待办（**DATA-1/2**）— 已按用户确认关闭
- 全文语言 polish 与过强词检查（**TEXT-1**）— 已完成
- CP2K relaxation 补实验（**EXP-1/2**）— 仍待做
- 第二节逐条跟踪表已细化（2026-08-28）

## 七、正文修改稿存档（已写入 `main.tex`，仅供回溯）

说明：以下 LaTeX 片段已全部写入 [`main.tex`](main.tex) 并编译通过。**日常修回跟踪以第二节为准**；本节仅保留原文备份，不再单独维护状态。

### 7.1 R1.2 — Methods / Computational Details（插在 Phase I 之前）

**插入位置**：[`main.tex`](main.tex) `\subsection*{Computational Details}` 首段之后、`In Phase I (Abstraction)` 之前。

**需同步**：[`reference.bib`](reference.bib) 新增 RCSR 条目 `okeeffe2008rcsr`：

```bibtex
@article{okeeffe2008rcsr,
  title={Reticular Chemistry Structure Resource (RCSR) database of nets},
  author={O'Keeffe, Michael and Peskov, Maxim A and Ramsden, Stuart J and Yaghi, Omar M},
  journal={Acta Crystallographica Section B: Structural Science},
  volume={64},
  number={3},
  pages={405--411},
  year={2008},
  publisher={International Union of Crystallography},
  doi={10.1107/S0108768108009956}
}
```

**正文段落**：

```latex
The workflow starts from a user-specified primitive quotient graph and one chosen net embedding together with degree-compatible building-unit pools. Benchmark topologies (\texttt{ffc}, \texttt{fjh}, \texttt{vba}, and related cases) were obtained from net embeddings associated with the Reticular Chemistry Structure Resource (RCSR)~\cite{okeeffe2008rcsr}. EPINET and ToposPro can supply equivalent inputs, but this study does not benchmark global net classification or cross-database equivalence. The method is not a net classifier: proper-rotation orbit reduction removes redundant building-block assignments within the fixed input, rather than deciding whether two differently named or differently embedded nets represent the same colored net.
% 中文翻译：该工作流程从用户指定的初基商图、一个选定的网络嵌入以及配位数相容的构建单元库出发。基准拓扑（\texttt{ffc}、\texttt{fjh}、\texttt{vba} 及相关算例）来自与网状化学结构资源（RCSR）~\cite{okeeffe2008rcsr} 关联的网络嵌入。EPINET 和 ToposPro 也可提供同类输入，但本文并不对全球网络分类或跨数据库等价性做基准测试。本方法不是 net classifier：真旋转轨道约化是在固定输入内去除冗余的 building-block 指派，而不是判断两个不同命名或不同嵌入的网络是否代表同一个 colored net。
```

**给您的直白解释**：审稿人说 “rotations 不能普遍解决 colored net 区分”——意思是 RCSR/EPINET 里两个 net 是否“同一个”不能单靠旋转解决。你们的数据从 RCSR 来，正好说明：**RCSR 负责提供已选 net；你们的方法负责在这个固定 net 上做 assignment 去冗余和 reconstruction**，不是替代 RCSR 做全网分类。

---

### 7.2 R3.1 — Introduction 对比表（插在 “present work addresses…” 段之后）

**插入位置**：Introduction 中 L85–86 “The present work addresses…” 段之后、“This interface is relevant to crystal engineering…” 之前。

**引导句 + 表**：

```latex
We do not replace ToBaCCo, MOFid/MOFkey, or downstream machine-learning screeners. Instead, the workflow reduces orientation-compatible redundancy before candidates enter those construction, identification, or screening pipelines. Table~\ref{tab:related_work} summarizes this positioning.
% 中文翻译：本文并不替代 ToBaCCo、MOFid/MOFkey 或下游机器学习筛选器。相反，这一工作流程是在候选进入这些构建、标识或筛选流程之前，先去除与取向相容的冗余。表~\ref{tab:related_work} 总结了这种定位关系。

\begin{table}[htbp]
\centering
\small
\caption{\textbf{Positioning relative to existing MOF generation and screening tools.} This work is scoped to candidate preparation on a fixed net rather than to replacement of downstream construction, identification, or screening workflows.}
% 中文翻译：与现有 MOF 生成和筛选工具的定位关系。本文聚焦于固定网络上的候选准备，而不是替代下游构建、标识或筛选流程。
\label{tab:related_work}
\setlength{\tabcolsep}{2.5mm}
\renewcommand{\arraystretch}{1.15}
\begin{tabular}{p{0.20\linewidth}p{0.24\linewidth}p{0.46\linewidth}}
\toprule
\textbf{Tool class} & \textbf{Representative examples} & \textbf{Relation to this work} \\
\midrule
Net-based assembly & ToBaCCo~\cite{colon2017topologically} & Builds MOFs from nets and building units; this work adds non-redundant ordered assignments and oriented reconstruction on a fixed net. \\
Identifier schemes & MOFid/MOFkey~\cite{bucior2019identification} & Supports search and deconstruction of existing or hypothetical MOFs; this work outputs traceable pre-relaxation starting models. \\
ML / free-energy screening & Refs.~\cite{niyongabo2025freeenergy,badrinarayanan2025mofgpt,Kim2024MOFFlowFM} & Downstream screening and generation; this work prepares a smaller, symmetry-aware candidate pool before those stages. \\
This work & Proper-rotation orbit enum.\ + Kabsch reconstruction & Removes reconstruction-compatible redundancy and records reconstruction diagnostics on fixed-net design spaces. \\
\bottomrule
\end{tabular}
\end{table}
```

---

### 7.3 R1.5 — Results 拓扑选取理由（插在 `\section*{Results}` 之后、第一个 subsection 之前）

```latex
The benchmark topologies were assigned distinct validation roles. \texttt{ffc} has a 7-node quotient graph and is used for exhaustive Burnside/CASM cross-checks and multi-composition DFT retention tests. \texttt{fjh} has a 28-node quotient graph and is used for two completed fixed-composition DFT comparisons across five $\rho$ strata. \texttt{vba} is used as a bounded-regime test when exhaustive enumeration is no longer practical.
% 中文翻译：各基准拓扑承担不同的验证作用。\texttt{ffc} 具有 7 节点商图，用于可穷举的 Burnside/CASM 交叉核对和多组成 DFT 保留效应检验。\texttt{fjh} 具有 28 节点商图，用于两个已完成的固定组成 DFT 比较，并覆盖 5 个 $\rho$ 层级。\texttt{vba} 用于可穷举枚举不再实际可行时的受限运行测试。
```

---

### 7.4 R3 minor 3 — Methods Building-Unit Inventory（插在 CASM 小节之前）

**同时修改 DFT 段**：将 “Because the frameworks combine transition metals with conjugated linkers” 改为 “Because the reported frameworks combine Cd or Cu nodes with conjugated organic linkers”，并更新对应 `%中文翻译`。

```latex
\subsection*{Building-Unit Inventory}
Table~\ref{tab:building_units} lists the node and edge building units used in the main-text enumeration, reconstruction, and DFT validations. Metal-containing nodes in the reported calculations are limited to Cd and Cu: Cd was used as the primary closed-shell node in the \texttt{ffc}/\texttt{fjh} DFT sets, while Cu appears in optional open-shell probe runs and in the \texttt{vba} visual-example pool. Organic nodes such as \texttt{3c\_B\_Ch} and \texttt{4c\_B\_Ch} are boron-based linkers rather than transition-metal SBUs.
% 中文翻译：表~\ref{tab:building_units} 列出了主文枚举、重构和 DFT 验证中实际使用的节点与边构建单元。报告的计算中，含金属节点仅限于 Cd 和 Cu：Cd 是 \texttt{ffc}/\texttt{fjh} DFT 集合中的主要闭壳节点，Cu 则出现在可选的开壳探针运行以及 \texttt{vba} 可视化示例池中。诸如 \texttt{3c\_B\_Ch} 和 \texttt{4c\_B\_Ch} 之类的有机节点是硼基连接体，而不是过渡金属 SBU。

\begin{table}[htbp]
\centering
\small
\caption{\textbf{Building units used in the main-text validation cases.} Edge linkers are treated as rigid connectors in reconstruction; metal nodes denote the transition-metal SBUs explicitly used in enumeration and/or DFT.}
% 中文翻译：主文验证算例中使用的构建单元。边连接体在重构中被视为刚性连接分子；金属节点表示在枚举和/或 DFT 中实际使用的过渡金属 SBU。
\label{tab:building_units}
\setlength{\tabcolsep}{2.5mm}
\renewcommand{\arraystretch}{1.15}
\begin{tabular}{llll}
\toprule
\textbf{Case} & \textbf{3-connected nodes} & \textbf{4-connected nodes} & \textbf{Edge linker} \\
\midrule
\texttt{ffc} enum./recon./DFT & \texttt{3c\_B\_Ch}, \texttt{3c\_triazine\_Ch} & \texttt{4c\_B\_Ch}, \texttt{4c\_Cd\_1\_Ch} & \texttt{1B\_4H\_Ch} \\
\texttt{fjh} enum./recon./DFT & \texttt{3c\_B\_Ch} & \texttt{4c\_B\_Ch}, \texttt{4c\_Cd\_1\_Ch} & \texttt{1B\_4H\_Ch} \\
\texttt{vba} visual example & --- & \texttt{4c\_B\_Ch}, \texttt{4c\_Cd\_1\_Ch}, \texttt{4c\_Cu\_1\_Ch} & \texttt{1B\_4H\_Ch} \\
\bottomrule
\end{tabular}
\end{table}
```

---

### 7.5 R1.2 + R3.3 — Limitations 追加段（插在 DFT 限制段之后、interpenetration 段之前）

```latex
A third limit concerns net-level equivalence and post-relaxation identity. Crystal rotations alone do not solve colored-net discrimination across databases or embeddings. That classification step must be completed before this workflow is applied. Distinct enumerated assignments may also relax to similar or identical geometries. Full cell and coordinate relaxation can merge nominally distinct candidates into the same basin. The single-point DFT scores reported here are not evidence that each symmetry class corresponds to a distinct relaxed material.
% 中文翻译：第三个限制涉及网络层面的等价性以及弛豫后的结构同一性。仅靠晶体旋转不能解决跨数据库或不同嵌入方式下的 colored net 区分问题；这一分类步骤必须在应用本文工作流之前完成。此外，名义上不同的枚举指派在弛豫后可能收敛到相似甚至相同的几何。完整的晶胞和坐标弛豫可能把不同候选并入同一 basin。本文报告的单点 DFT 分数不能作为每个对称性类别都对应互异弛豫后材料的证据。
```

---

### 7.6 Data availability 措辞（公开链接）

当前已使用 definitive public URL：

```latex
The code and supporting data required to reproduce the findings of this study are publicly available at GitHub: https://github.com/Buct454/mof-enumeration-reconstruction.git. The package contains ...
% 中文翻译：用于复现本文结果的代码和支撑数据已公开存放于 GitHub：https://github.com/Buct454/mof-enumeration-reconstruction.git。该数据包包含 ...
```

---

### 7.7 写入后编译检查

```bash
cd paper && latexmk -pdf main.tex
```

确认 `tab:related_work`、`tab:building_units`、`okeeffe2008rcsr` 引用无 undefined。

---

## 八、待实验 / 待 push TODO

| TODO ID | 内容                                                                 | 关联意见         | 状态 | 备注                                                                            |
| ------- | -------------------------------------------------------------------- | ---------------- | ---- | ------------------------------------------------------------------------------- |
| EXP-1   | 少量 CP2K `CELL_OPT`：Fig.2 三个代表 CIF + 2–3 个 DFT minima      | R3.3             | 待做 | 脚本：`dft_validation/run_relaxation_validation.py`；先续跑 `ffc_comp_a_11` |
| EXP-2   | 仅**收敛且 topology-retained** 的结果写 SI 小表                | R3.3             | 待做 | 未收敛不能写进正文                                                              |
| MTV-1   | **结构选取**：从现有 CO2_N2 `batch_002_clean_pool_5000` 中选 3 个 matched `ffc` comparison groups；每组固定 composition、\(\rho\)、symmetry stratum，并含 4 个 ordering | R1.5, R3.5       | 已完成 | 共 12 个候选级标签；主文表格按 group 汇总，完整 ID 在 CSV 中 |
| MTV-2   | **GCMC 前结构准备**：复用 batch_002 clean-pool 中已完成的 Zeo++ 描述符和结构筛选字段 | R1.5, R3.5       | 已完成 | 主文报告 PLD/LCD；本地数据锚点：CO2_N2 `data/interim/co2n2_ffc_screening/batch_002_clean_pool_5000/` |
| MTV-3   | **性质计算（首选 GCMC）**：复用已有 CO2/N2 GCMC-derived labels，统一 298 K、1.0/0.1 bar、15/85 mixture、RASPA2 设置 | R1.5, R3.5       | 已完成 | 主指标：CO2 working capacity、CO2/N2 selectivity、regenerability、APS；不写 stability/synthesizability |
| MTV-4   | **SI/Data package 整理**：主文 3 行汇总表 + compact CSV 数据锚点 | R1.5, R3.5       | 已完成 | 新增 `data_availability/mof-enumeration-reconstruction/tables/ffc_co2n2_ordering_check.csv`，并更新 manifest / figure-table mapping；CSV 含 12 行完整候选数据 |
| MTV-5   | **正文/Rebuttal 回写**：Results 增 fixed-composition ordering check；Discussion/Conclusion/Methods/Data availability 同步 | R1.5, R3.5       | 已完成 | 仅写 screening-stage / adsorption-derived screening labels；不夸大 sequence-dependent synthesis |
| DATA-1  | 推送 `data_availability/mof-enumeration-reconstruction/` 到 GitHub | R1.1, R3.2, R3.4 | 已关闭 | 用户已确认数据包已 push 到公开仓库；正文 Data availability 已改为 definitive public URL；MTV-4 如产生新文件，可后续补入 |
| DATA-2  | push 后核对链接可访问、manifest/CSV 与主文表一致                     | R1.1, R3.2       | 已关闭 | 用户已确认公开仓库状态；主文当前使用公开 GitHub URL |
| SI-1    | 从 representative CIF 导出 SI 全尺寸结构图                           | R3 minor 2       | 可选 | Fig.2 主文已放大                                                                |
| AUDIT-1 | `vba` 完整 reconstruction audit（若需与 ffc/fjh 同表）             | R3.4             | 可选 | 当前 vba 仅 1/10 完整通过                                                       |
| TEXT-1  | 全文语言 polish（Abstract/Discussion 压缩长句；过强词检查）          | R1.4, R1.6       | 已完成 | 已完成一轮全稿中性化、降调和中文翻译同步；最终 rebuttal 前可再做只读复扫 |
| REBUTTAL-1 | 逐点回复信（nature-response）                                     | 全部             | 待 MTV-5 / DATA 定稿后 | 可引用 §7 rebuttal 口径                                                         |

**方案 C 执行记录（MTV-*）**

#分优先级,重要问题处理，重点，贡献不超过三条，收紧，可复现性*


1. **MTV-1 选结构**
   - 拓扑：`ffc`
   - 来源：CO2_N2 `batch_002_clean_pool_5000`
   - 固定组成组：3 个 matched groups，分别为 A/B/C。
   - 数量：每组 4 个 symmetry-distinct ordering，共 12 个候选级标签。
   - 候选：A = `00001/00004/00008/00011`；B = `00002/00005/00007/00010`；C = `00013/00016/00019/00022`。

2. **MTV-2 结构准备**
   - 复用 existing clean-pool 筛选表；主文只报告 PLD/LCD，避免把孔结构描述符解释成最终吸附性能。

3. **MTV-3 性质计算**
   - 复用已有 RASPA2/GCMC-derived labels：298 K，CO2/N2 = 15/85，adsorption/desorption pressure = 1.0/0.1 bar。
   - 成功结果：3 个 matched groups 的 APS 组内范围分别为 0.595–0.817、0.588–0.760、0.626–0.839；CO2 working capacity 组内范围分别为 0.345–0.411、0.341–0.386、0.349–0.422 mmol/g。

4. **MTV-4 → MTV-5 写回**
   - 正文：在 Results DFT 图后新增 `Matched Ordering Check for CO2/N2 Screening` 小节和表 `tab:ffc_co2n2_ordering`。
   - Methods：新增 CO2/N2 screening labels 小节，定义 selectivity、working capacity、regenerability 和 APS。
   - Data package：新增 compact CSV，并更新 manifest 与 figure/table mapping。
   - Rebuttal R1.5/R3.5：引用该主文表作为 crystal-engineering 接口的最小实证。

**MTV-1 结构候选来源（已确认）**

- 从 CO2_N2 `batch_002_clean_pool_5000` 的 `candidates.csv`、`descriptors.csv` 和 `gcmc_labels.csv` 合并筛选；
- 选中 3 个 matched comparison groups；每组内部 composition、\(\rho\)、symmetry stratum 固定，且 4 个 ordering 的 GCMC 均 succeeded；
- 该选择避免把观察结果解释为 composition、\(\rho\) 或 symmetry-stratum 变化造成的差异。

**与现有代码的关系**

- 枚举/选结构：`enumerator.py`, `dft_validation/select_structures.py`, `dft_validation/config.py`
- 重构/验证：`assembler.py`, `experiment/geometric_reconstruction_summary.py`, `structure_validator.py`
- 吸附候选导出（可参考，非 MTV 专用）：`adsorption_screening/export_ffc_co2n2_candidates.py`, `export_ffc_co2n2_clean_pool.py`
- GCMC：本次 MTV 小实验复用外部 CO2_N2 项目中已完成的 RASPA2 标签；主文 Methods 已记录输入、力场和条件

**EXP-1 建议批次**（来自 `codex/NEXT_ACTIONS.md`）：

- Fig.2 代表：`ffc_comp_b_21`, `fjh_comp_b_1`, `vba_rho1`
- DFT minima：`ffc_comp_a_11`, `ffc_comp_b_21`, `fjh_comp_a_12`
- 预设：`--preset reviewer7` 或 `--only` 单结构

**Rebuttal 口径（R1.2）**：

> We agree that crystal rotations alone do not universally solve colored-net discrimination across databases. Our method is scoped to enumeration and reconstruction on a fixed quotient graph obtained from RCSR net embeddings. We do not claim to replace RCSR/EPINET/ToposPro net analysis.
