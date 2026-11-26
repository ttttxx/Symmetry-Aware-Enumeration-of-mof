# main_MATCH.tex 联动改稿拟稿（2026-11-16）

状态：供逐条核对的拟替换文字；尚未写入主稿，不代表作者已经接受以下表格布局或英文措辞。原主稿以本次读取的文件为准。

本轮范围：普通枚举算法、相关正确性论证、实现与复杂度说明、同群性能比较，以及它们牵涉的摘要/贡献/图注/模板库/结构审计/配额附录/结论。结构构建方法、吸附和 DFT 数值保持现有证据范围。

## 2026-11-23 八处统一方案已批准并落稿（最新状态）

作者“OK 改”批准了紧前一版取消“普通调度”名称的八处联动方案。现已修改main_MATCH.tex的任务组织、完整覆盖说明、配额可选处理、同层覆盖及受限条件、Algorithm 1、配额附录证明，并新增实现附录短段。以完整枚举为主流程，配额为大空间下限制输出的补充；主文不介绍Trie布局。之前含ordinary schedule的TeX候选已由本次主稿取代，不再作为待应用版本。实验结果与表格、摘要、结论及数据声明本轮保持不变；第10项等其余讨论继续待定。作用范围核验见MATCH_UNIFIED_APPLIED_SCOPE_2026-11-23.json，修改前源码快照见snapshots/match_before_unified_method_20261123.tex。


## 2026-11-23 最新定位与审批（优先于以下历史方案）

完整枚举是主流程，精确ρ配额是规模过大时限制输出的补充。不为主流程另起“普通调度/ordinary schedule”名称，不按两套并列方法展开。上一轮完整方案第2项已批准并落稿；第4–8项内容认可，但涉及上述名称或主次定位的位置先待定，待统一文字后再修改。此要求不代表启用自动超限切换配额，代码启用条件保持现状。


## 2026-11-23 本轮更新：第1/14项已批准，统一方法与第10项继续审阅

- 已将上一轮第1项两句摘要和第14项结论首段写入main_MATCH.tex，中文同步。引言贡献未另改；方法段、Algorithm 1、配额证明、CASM主表及vba表保持本轮开始时状态。完整编译32页，六处原有overfull不变。
- 第4、6–8项采用同一优化DFS核心、普通/精确ρ两种调度的方向获认可；完整文字仍供本轮审阅，尚未落稿。完整2.2与统一Algorithm 1位于 `MATCH_UNIFIED_METHOD_PROPOSAL_2026-11-23.tex`；联动配额证明位于 `MATCH_QUOTA_PROOF_PROPOSAL_2026-11-23.tex`；移到实现附录的短段位于 `MATCH_IMPLEMENTATION_PROPOSAL_2026-11-23.tex`。主文不新增Trie结构细节。固定点域/规范拒绝公式保留，配额过滤先于完整键保留；覆盖证书按接收域限定，受限条件集中一处。
- 第10项继续讨论：正文使用enumeration time，计时范围仅在表注/附录定义一次。记录一致用于验证正确性；同群、同输入和同输出字段是预先统一的比较条件，二者不写成因果推导。性能结果用于说明被测较大空间中完整候选库生成的实际效率；子群调度的分层输出用途由vba实验支撑。
- 撤回“固定成本超过剪枝收益”的具体因果句。已有计时未完整分解全部成本及逐项消融，不把小例差距全归于某一项。讨论保留子群任务的理由：它们组织固定点搜索，并支持精确ρ配额下的分层输出；CASM计时与vba层分布分别评价效率及输出控制。

### 第10项待讨论正文

Both workflows used the same group Γ and degree-compatible label pools and retained complete representatives, stabilizer orders, and orbit sizes. Independent orbit-count and coverage checks verified the returned candidate sets.

译文：两种流程采用相同的Γ群及度数相容标签池，并保留完整代表、稳定子阶数和轨道大小。独立的轨道计数与覆盖核查验证了返回的候选集合。

For the tested larger assignment spaces, the proposed implementation generated the complete representative set in less time than the matched CASM workflow. These results support its practical use for exhaustive candidate generation at the tested scales. The subgroup tasks also organize fixed-point searches for exact-ρ quotas, whose output across symmetry strata is examined in the bounded-vba comparison.

译文：对于被测的较大指派空间，本文实现生成完整代表集合所需时间少于采用匹配条件的CASM流程。这些结果支持其在被测规模下用于实际的穷尽候选生成。子群任务还为精确ρ配额组织固定点搜索，其输出在各对称层间的分布由受限vba对照检验。

CASM advances a counter over complete assignments. The subgroup schedule used here adds task-construction work and provides fixed-point domains for organizing the search by symmetry. This organization supports both exhaustive enumeration and bounded output with prescribed per-stratum quotas.

译文：CASM通过计数器逐个生成完整指派。本文采用的子群调度增加了任务构造工作，同时提供固定点域以按对称性组织搜索。这种组织方式支持穷尽枚举及具有指定分层配额的受限输出。

说明：第三段描述结构成本及采用此组织方式的目的，未声称已完成小例耗时差距的完整归因。最终结果正文是否保留这一段，继续由作者评议。

---

## 2026-11-23 逐项处理状态（本节优先于下方历史拟稿）

- 作者批准术语；论文仅落实第2项（删除末句受限说明）、3、5、9、11、13项。第3项的叶节点/游标同步到原Algorithm 1，未套用第6项整套重写。已批准正文采用既有global membership check表述，没有扩写正文Trie实现。
- 第1、4、10、12、14项继续讨论；第7、8、15项在配额代码完成后讨论，未据此自动改稿。摘要/贡献、CASM主表、vba表、结论和数据声明保持修改前文字。
- 第6项代码已经实现：精确ρ配额共享优化Python/C++搜索核，按目标稳定子阶数过滤，覆盖证书按层隔离；ABI 4。211项测试通过，论文vba实际池18次独立测量的完整有序输出、各层计数与组成向量均一致。此前拟稿中“新后端未整合配额”“精确ρ配额独立旧Python执行”等表述已被代码进展取代，不能继续作为拟落稿文字。
- 新增CASM附录使用冻结ABI-3数据，当前ABI-4配额新时间仅保存在验证记录中。主表未增加速度比列，附表也不增加该列。

### 继续讨论的替换文字（尚未应用）

**第1项摘要概括：**

We organize fixed-point searches through the subgroup lattice and use depth-first search with canonical prefix comparisons to identify orbit representatives. Completed search regions are reused across tasks to avoid repeated exploration.

译文：我们通过子群格组织固定点搜索，并采用带规范前缀比较的深度优先搜索识别轨道代表。已完成的搜索区域在任务间复用，以避免重复探索。

理由：概括搜索组织、规范判据和跨任务复用；完整键索引属于实现说明，摘要不突出Trie。

**第4项正文：** 不增加原第4项的独立Trie段，沿用规范代表定义及已批准完备性段的global membership check。数据结构和数组布局保留在代码/复现材料。

**第10项CASM结果解释：**

Under the same group action, the two workflows returned identical representative records, including stabilizer orders and orbit sizes. The timing results therefore assess the cost of generating the same symmetry-distinct candidate library. For the tested larger assignment spaces, the lower wall times support exhaustive enumeration with the proposed search procedure. The small-case timings are consistent with the fixed costs of subgroup-task preparation and search setup outweighing the savings from prefix pruning.

译文：在相同群作用下，两种流程返回了相同的代表记录，包括稳定子阶数和轨道大小。因此，计时结果评价的是生成同一个对称互异候选库的成本。对于被测的较大指派空间，更低的壁钟耗时支持采用所提出的搜索过程进行穷尽枚举。小算例的计时与以下解释相符：子群任务准备及搜索初始化的固定成本超过了前缀剪枝节省的工作。

说明：最后一句属于基于阶段计时的解释，尚未通过完整准备开销消融量化；措辞保留consistent with。主表只列双方时间，不单列速度比。

**第12项：** 原优化前后vba时间叙述撤回，不加入结果正文。配额实现已完成，待讨论后统一使用新配额与匹配普通DFS的实验协议。代码回归报告可保留优化前后时间，论文正文聚焦各对称层的候选覆盖。

**第14项结论候选：**

In exhaustive runs, the method returns one representative per Γ-orbit. The matched-group benchmarks support its use for generating symmetry-distinct candidate libraries over the tested assignment spaces. Selected representatives are then used to construct pre-relaxation periodic models for subsequent screening.

译文：在穷尽运行中，该方法为每个Γ轨道返回一个代表。同群基准支持将其用于被测指派空间中对称互异候选库的生成。随后，以选定代表构建弛豫前周期模型，供后续筛选使用。

理由：结论保留方法产出、实用性及模型用途；小例的成本解释集中在CASM实验部分。

**配额实现后第6–8项的联动方向：** 普通调度与精确ρ配额采用同一优化DFS核心；Algorithm 1可以保留可选配额分支。需要一并核对2.2的任务顺序/精确层过滤、同层覆盖证书、Algorithm 1、配额附录证明及vba方法说明。截止条件集中交代，避免在各段重复。具体正文等作者逐项认可后再落稿。

---

## 术语与符号

| 统一英文 | 中文 | 用法 |
|---|---|---|
| node assignment / coloring c | 节点指派 / 着色 c | 沿用前文定义；不称原子构型 |
| canonical representative | 规范代表 | Γ 轨道的字典序最小元 |
| canonical prefix comparison / rejection | 规范前缀比较 / 拒绝 | 群像字典序比较，与覆盖缓存查询分开 |
| completed-task coverage | 已完成任务覆盖 | 来源任务整个固定点域的规范搜索已结束 |
| completed-prefix cache | 已完成前缀缓存 | 来源任务已完成子树的覆盖证书 |
| representative Trie T_R | 代表 Trie T_R | 完整规范键的精确成员索引 |
| active comparison set A | 活跃比较集 A | 保留尚未确定严格大于当前指派的非恒等操作 |
| exact-ρ quota schedule | 精确 ρ 配额调度 | 独立任务次序和精确稳定子过滤，不声称已经使用新 C++ 后端 |
| enumeration wall time | 枚举壁钟耗时 | 任务准备至结果排序/解码；不称全部 C++ 核时间 |

以下暂用小节标签 `app:enum_implementation`、`app:enum_benchmark`、`alg:quota`，落稿时核对编号和交叉引用。正文现有 `alg:enum`、`tab:casm`、`app:quota_proof` 保留。

## 1. 摘要枚举句、第一项贡献和流程图(b)

**摘要替换一句：**

We organize fixed-point searches through the subgroup lattice and combine canonical prefix comparisons with a global representative Trie to return distinct orbit representatives.

**译文：** 我们通过子群格组织固定点搜索，并结合规范前缀比较与全局代表 Trie，返回互不重复的轨道代表。

**第一项贡献替换为：**

We organize deterministic fixed-point searches through the subgroup lattice and exclude regions already covered by completed tasks. Canonical prefix comparisons and a global representative Trie yield one representative per orbit under the chosen symmetry group when the search is exhaustive.

**译文：** 我们通过子群格组织确定性的固定点搜索，并排除已完成任务覆盖的区域。在穷尽搜索条件下，规范前缀比较与全局代表 Trie 为所选对称群下的每个轨道返回一个代表。

**流程图 caption 的(b)替换为：**

(b) subgroup fixed-point tasks return canonical representatives using the ordinary schedule or the exact-ρ quota schedule;

**译文：** (b) 子群固定点任务采用普通调度或精确 ρ 配额调度，返回规范代表；

**理由：** 保留作者确定的 subgroup lattice；明确搜索中判规范与全局查重的分工。完整覆盖条件继续在贡献项中保留；图注的两种调度与后面的实现分工一致。摘要不增加单个最大速度比。

## 2. 2.2 子群任务段后增加已完成区域排除

Different subgroup tasks can have overlapping fixed-point domains. After a task has been searched to completion, its vertex-orbit partition is recorded as a set of equal-label constraints. A later branch is rejected when every completion satisfies the constraints of an already completed task. This coverage test can discard an entire branch before a complete assignment is formed. A task stopped by a time or output limit contributes only the subtrees whose searches have finished.

**译文：** 不同子群任务的固定点域可以相交。一个任务搜索完成后，其顶点轨道划分被记录为一组标签相等约束。若后续分支的每个补全均满足某个已完成任务的约束，则拒绝该分支。这项覆盖检查可以在形成完整指派之前排除整个分支。因时间或输出限制停止的任务，仅提供其中已经完成搜索的子树。

**理由：** 原文只说明相同划分合并，没有交代不同划分间的交叠如何处理。保留现有任务顺序段，不把 Fix(H) 说成天然互斥。最后一句概括可选前缀缓存，适用条件放第4项与附录。

## 3. 规范前缀判据后增加比较游标和叶节点处理

For each operation in the active comparison set, the search stores the end of the vertex prefix already confirmed to be equal. Descendant states resume comparison from that position, and backtracking restores the comparison state. At a complete surviving leaf, the coloring is already the canonical representative. The active set A then contains exactly the nonidentity operations that fix this coloring, giving

s(c) = |Stab_Γ(c)| = 1 + |A|,     |Γ·c| = |Γ| / s(c).

The leaf and these two quantities are passed directly to the representative Trie.

**译文：** 对活跃比较集中的每个操作，搜索保存已经确认相等的顶点前缀的结束位置。后代状态从该位置继续比较，回溯时恢复比较状态。在完整且未被拒绝的叶节点，当前着色已经是规范代表。此时，活跃集 A 恰好由保持该着色不变的非恒等操作组成，因此 s(c)=|Stab_Γ(c)|=1+|A|，轨道大小为 |Γ|/s(c)。完整叶及这两个数值直接交给代表 Trie 处理。

**理由：** 与现有规范拒绝公式衔接；说明为什么能够复用稳定子信息。Algorithm 1 中再次求 κ(c) 和重新扫描全群的步骤必须同步删除。保留 κ 的数学定义，省去的是执行时重复计算。

## 4. 完整代表 Trie 与可选前缀缓存

**正文增加：**

The global representative Trie T_R indexes complete canonical assignments. Each key is compared exactly; a new key adds one output record containing the assignment, stabilizer order, and orbit size. The implementation uses a path-compressed byte Trie stored in contiguous node arrays, with edge labels referencing the retained assignment keys. A separate completed-prefix cache records finished subtrees of interrupted tasks. A cached prefix can exclude a later branch only when every completion obeys the source task's equal-label constraints and uses labels admitted by that task. Implementation details are given in Supplementary Section app:enum_implementation.

**译文：** 全局代表 Trie T_R 对完整规范指派建立索引。每个键均进行精确比较；新键对应一条输出记录，包含指派、稳定子阶数和轨道大小。实现采用连续节点数组存储路径压缩的字节 Trie，边标签引用已经保留的指派键。另设的已完成前缀缓存记录中断任务中已经完成的子树。仅当后续分支的每个补全都满足来源任务的标签相等约束，且所用标签均属于来源任务允许的范围时，缓存前缀才可用于排除该分支。实现细节见补充材料 app:enum_implementation。

**附录增加：**

A completed-prefix entry stores the source partition, admissible label domains, and a variable prefix whose canonical search has finished. The open path at interruption remains searchable. A later task uses an entry only after its remaining assignment domain is contained in the source task's domains and its completions are proved to satisfy the source equalities. Cache capacity limits the number of retained certificates. Once capacity is reached, further certificates may be omitted and enumeration continues with the remaining pruning tests.

**译文：** 每条已完成前缀记录保存来源划分、允许标签域，以及规范搜索已经完成的变量前缀。中断时的开放路径仍可继续搜索。只有当后续任务的剩余指派域包含于来源任务的允许域，且能够证明其补全满足来源等式时，才使用该记录。缓存容量限制保留的证书数量；容量用尽后可以不再保存额外证书，枚举仍使用其余剪枝检查继续进行。

**理由：** 完整键查重与搜索区域覆盖具有不同含义，需要两套状态。附录交代安全查询条件和容量边界；不宣称该缓存已在本轮基准贡献提速。

## 5. 替换 Completeness 全段

The identity task has domain Fix({e})=X, so every canonical representative lies in its search space. Canonical prefix rejection removes only branches whose completions have a smaller Γ-image. Coverage rejection removes a branch only after its canonical search has been completed in an earlier task or subtree. Each canonical representative is therefore either retained earlier or remains reachable in the identity task. Exact insertion into T_R records it once. With output limits and task timeouts disabled and all tasks searched to completion, the algorithm returns every member of Ω exactly once.

**译文：** 恒等任务的搜索域为 Fix({e})=X，因此每个规范代表均位于其搜索空间中。规范前缀拒绝仅移除补全存在更小 Γ 群像的分支。覆盖拒绝仅在该分支的规范搜索已由先前任务或子树完成后才排除它。因此，每个规范代表要么已在此前保留，要么仍能在恒等任务中到达。向 T_R 的精确插入确保它只记录一次。关闭输出限制与任务超时，并使全部任务搜索完成时，算法恰好返回 Ω 的每个成员一次。

**理由：** 原文依赖叶节点再规范化的表述被替换；增加新覆盖剪枝的正确性环节。普通流程无精确层过滤，配额正确性继续由独立附录论证。

## 6. Algorithm 1 整体替换：普通调度

**标题：** Subgroup-Lattice Enumeration with Canonical Prefix Comparisons and Coverage Reuse

**译文：** 结合规范前缀比较与覆盖复用的子群格枚举。

下面是完整逻辑拟稿；排版时将主流程放正文，DFS 或覆盖证书辅助细节可移附录，正文 Algorithm 1 仍明确引用。避免让长算法再次成为过高的浮动体。

```text
Require Γ, degree-compatible pools {C_d}, fixed vertex and label orders,
        optional output/time limits, optional completed-prefix caching
Ensure  distinct canonical representatives with stabilizer orders and orbit sizes

Encode labels as integers preserving their comparison order
R ← empty output records; T_R ← empty complete-key representative Trie
C ← empty completed-task coverage; F ← empty completed-prefix cache
Q ← subgroup-orbit tasks including Γ and {e}
Merge identical partitions as in the ordinary schedule
Order tasks by increasing number of vertex orbits, retaining tie order

for each partition P in Q:
    if the global output limit is reached: break
    Compile completed-task equalities from C for P
    Prepare admissible completed-prefix queries from F
    Initialize c_p ← unassigned, A ← Γ \ {e}, equal-prefix cursors ← 0
    Initialize task limits and completed-sibling progress
    state ← CanonicalDFS(0, c_p, A, cursors, P)
    if state = CLOSED and no task cutoff was triggered:
        Register P in C
    else if completed-prefix caching is enabled:
        Record only finished sibling subtrees from the interrupted frontier in F

Stable-sort R by decreasing stabilizer order; decode labels; return R

CanonicalDFS(k, c_p, A, cursors, P):
    if a time or output limit is reached:
        Save the interrupted frontier; return OPEN
    if completed-task coverage contains every completion: return CLOSED
    if a compatible completed prefix contains every completion: return CLOSED
    if k = |P|:
        s ← 1 + |A|
        if T_R.insert(c_p) reports a new complete key:
            Append (c_p, s, |Γ|/s) to R
        return CLOSED

    for each label a in the supplied pool order for block P[k]:
        if a time or output limit is reached:
            Save the interrupted frontier; return OPEN
        Assign a to all vertices in P[k]
        Copy the branch comparison state
        Advance each active operation from its equal-prefix cursor
        Stop each comparison at an unassigned pair or its first unequal pair
        If an image is smaller: mark the branch canonically rejected
        If an image is larger: remove its operation from the child's active set
        Otherwise retain the operation and its cursor
        child_state ← CLOSED if canonically rejected,
                      otherwise CanonicalDFS(k+1, child state, P)
        Restore assignments and comparison state
        if child_state = OPEN: return OPEN
        Mark this label choice as finished
    return CLOSED
```

**中文逐步对应：** 保序整数编码；建立结果表、完整键 Trie、完成任务覆盖表和可选前缀缓存；按原普通次序调度子群任务；将先前覆盖条件编译到当前任务；初始化活跃集与游标；DFS 先查限制及覆盖，再赋标签并推进规范比较；完整叶直接以 1+|A| 求稳定子并精确插入；完整任务登记整域，中断任务只记录已完成兄弟子树；最后按稳定子阶数稳定排序、解码并返回。

CLOSED 表示当前任务约束下的规范搜索已完成或被安全排除；OPEN 表示因限制停止、仍有未处理路径。

**从原 Algorithm 1 联动移出的内容：** 精确层配额 B_s、目标阶 s、层计数 N_s、按阶数递减调度及精确阶数过滤，移入附录 Algorithm `alg:quota`。原普通调度的同分顺序和标签池遍历顺序保留；比较标签的全局次序与实际标签遍历顺序区分。

## 7. 精确 ρ 配额的正文过渡、算法与附录边界

**正文配额段增加/调整为：**

For bounded output across symmetry strata, a separate exact-ρ quota schedule assigns a budget B_s to each Ω_s. It orders tasks by decreasing subgroup order and retains a representative from an order-s task only when its full stabilizer order is s. The quota implementation uses canonical prefix rejection and complete-representative deduplication. A completed search for one stabilizer stratum supplies no certificate that the entire fixed-point domain has been covered. Its task construction and coverage guarantee are given in Supplementary Section app:quota_proof and Algorithm alg:quota.

**译文：** 为在受限输出中保留各对称层的候选，独立的精确 ρ 配额调度为每个 Ω_s 分配预算 B_s。任务按子群阶数递减排序，阶数为 s 的任务仅保留完整稳定子阶数恰为 s 的代表。配额实现采用规范前缀拒绝和完整代表去重。对某个稳定子层完成搜索，并不能作为整个固定点域已覆盖的证书。任务构造和覆盖保证见补充材料 app:quota_proof 及算法 alg:quota。

**附录首段增加：**

Algorithm alg:quota specifies the exact-ρ quota schedule used in the bounded-vba comparison. It retains the Python canonicalization and leaf deduplication procedure and applies the exact-order acceptance rule. Completed-task and completed-prefix coverage reuse in Algorithm alg:enum applies to the ordinary schedule described there.

**译文：** 算法 alg:quota 给出受限 vba 对照中采用的精确 ρ 配额调度。该实现保留 Python 规范化及叶节点去重过程，并施加精确阶数接收规则。算法 alg:enum 中的已完成任务与已完成前缀覆盖复用适用于该算法所述的普通调度。

**拟新增配额算法（保留原实现规则）：**

```text
Require Γ, degree-compatible pools, per-stratum budgets B_s,
        task timeout, total output cap
Build all subgroup partitions, deduplicated within each subgroup order
R ← empty records; K_seen ← empty complete-canonical cache
for s in decreasing subgroup order:
    N_s ← 0
    for each order-s partition P in lexicographic partition order:
        if N_s = B_s or global cap is reached: break
        Run canonical-prefix DFS over Fix(P) with the task deadline
        At each surviving complete leaf c:
            c* ← κ(c)
            if c* is already in K_seen: skip this leaf
            Insert c* into K_seen
            t ← |Stab_Γ(c*)|
            if t = s:
                Append (c*, t, |Γ|/t) to R; N_s ← N_s + 1
            Stop this task/layer when its corresponding limit is reached
Sort retained records by stabilizer order; return R
```

**译文：** 按子群阶数递减组织任务，在同阶内合并相同划分并按字典序搜索；叶节点沿用规范化和全局完整键查重，仅接收完整稳定子阶数等于任务阶数的代表；仅被接收的代表消耗该层预算，并遵循任务超时和总输出上限。

**保留项：** Ω_s、ρ=s/|Γ|、前 min{B_s,|Ω_s|} 个代表的条件保证，以及附录原证明/公式均保留。原证明中的叶节点 κ(c) 用语对应配额实现，不能随普通算法一起全局删去。

## 8. 实现及时间/空间成本

**正文增加：**

Python constructs the subgroup tasks and prepares their coverage constraints. A C++ implementation executes each task's DFS, canonical comparisons, stabilizer extraction, and representative-Trie insertion. Results are transferred in batches and decoded and sorted in Python. Labels are encoded in a common integer order while each task retains its prescribed label traversal order. The exact-ρ quota runs use the separate Python implementation described in Supplementary Section app:quota_proof.

**译文：** Python 构造子群任务并准备覆盖约束。C++ 实现执行各任务的 DFS、规范比较、稳定子提取和代表 Trie 插入。结果分批传回，再由 Python 解码和排序。标签按统一次序编码为整数，各任务保留规定的标签遍历顺序。精确 ρ 配额运行采用补充材料 app:quota_proof 所述的独立 Python 实现。

**附录成本段：**

For n vertices and r retained representatives, the path-compressed representative Trie has O(r) nodes, while the stored assignment keys require O(nr) space. Each nonroot internal node branches into at least two children, giving at most 2r+1 Trie nodes. Total memory also includes output records, subgroup partitions, comparison states, coverage constraints, and transfer buffers. Enumeration remains exponential in the worst case; the identity task has the full assignment domain X. The optimizations reduce repeated work and storage overhead within this search.

**译文：** 对 n 个顶点和 r 个保留代表，路径压缩代表 Trie 有 O(r) 个节点，而存储指派键需要 O(nr) 空间。每个非根内部节点至少分出两个孩子，因此 Trie 节点数不超过 2r+1。总内存还包括输出记录、子群划分、比较状态、覆盖约束和传输缓冲。枚举在最坏情况下仍为指数级；恒等任务具有完整指派域 X。这些优化减少该搜索中的重复工作和存储开销。

**实现附录还须明确的配置：** cpp 后端、native_store=trie、prefix_trie=True、ABI 3；普通穷尽测试关闭任务超时/总数上限/恒等任务新增配额。配额测试配置仍对应原 vba 表。前缀缓存默认容量50,000；本轮所有正式测量 queries/hits均为0。宏观数学正文无需写 ABI 或命令行参数。

## 9. Results 开头、标签池及构建审计来源

**Results 首段替换为：**

The ffc and fjh nets were used for exhaustive enumeration benchmarks, and vba was used to examine bounded output across symmetry strata. The timing benchmarks varied the sizes of the degree-compatible node-label pools; their exact membership is recorded in Supplementary Section app:enum_benchmark. Models for construction audits and property screening used the separately specified building-unit pools in Supplementary Table tab:building_units, with 1B_4H_Ch fixed as the linker on every edge.

**译文：** ffc 和 fjh 网络用于穷尽枚举基准，vba 用于考察受限输出在各对称层间的分布。计时基准改变度数相容节点标签池的大小，其具体成员记录于补充材料 app:enum_benchmark。用于构建审计和性质筛选的模型采用补充表 tab:building_units 中另行指定的构建单元池，每条边均固定使用 1B_4H_Ch 连接体。

**构建审计开头补充/替换为：**

For the construction audit, the ffc pool comprised 3c_B_Ch and 3c_triazine_Ch at three-connected vertices, and 4c_B_Ch and 4c_Cd_1_Ch at four-connected vertices. The fjh pool used 3c_B_Ch at three-connected vertices and the same two four-connected templates. All 42 ffc representatives and a fixed-seed sample of 30 out of 394 fjh representatives were assembled successfully.

**译文：** 构建审计的 ffc 池在三连接顶点使用 3c_B_Ch 和 3c_triazine_Ch，在四连接顶点使用 4c_B_Ch 和 4c_Cd_1_Ch。fjh 池在三连接顶点使用 3c_B_Ch，四连接顶点使用相同的两种模板。ffc 的全部42个代表，以及从 fjh 的394个代表中按固定种子抽取的30个代表均成功组装。

**附录审计协议同步：** 增加“采用本节明确列出的审计池”，或直接交叉引用本段，保留抽样种子/层化规则及原72个模型。正文紧接此段的 CIF、topology recovery、short-contact 等已有结论保留。

**理由：** 新计时池按文件名排序取前若干标签，不是构建审计的同一化学模板集合。相同基数的池可通过标签重命名给出相同轨道计数；不能将新池下同样的42/394误写成已经审计了新模板。这里是本轮必须联动的样本来源修正。

## 10. CASM 正文、表2与性能结论

**小节标题建议：** Exact Enumeration Checks and Matched-Group Performance

**译文：** 精确枚举核查与同群性能比较。

**验证段替换为：**

For each exhaustive case, the representative count agreed with an independent Burnside calculation, and the orbit sizes summed to |X|. CASM 2.3.1 was evaluated on the same quotient graph and node-label pools using a subgroup whose site permutations matched Γ. After conversion to a common canonical convention, the complete representative records agreed, including stabilizer orders and orbit sizes.

**译文：** 对每个穷尽算例，代表数均与独立 Burnside 计算一致，轨道大小之和等于 |X|。CASM 2.3.1 在相同商图和节点标签池上运行，其指定子群的位点置换与 Γ 一致。转换到统一的规范代表约定后，全部代表记录一致，包括稳定子阶数和轨道大小。

**主表 caption：**

Enumeration under the same permutation group Γ. Both methods retain complete representatives, stabilizer orders, and orbit sizes. Wall times are medians of three independent runs and include task preparation, enumeration, result storage, and sorting. Graph and group construction, dependency imports, and verification are excluded. Speedup is the CASM wall time divided by the wall time of this work.

**译文：** 相同置换群 Γ 下的枚举。两种方法均保留完整代表、稳定子阶数和轨道大小。壁钟耗时为三次独立运行的中位数，包括任务准备、枚举、结果保存和排序。图与群构建、依赖导入和验证不计入。速度比为 CASM 壁钟耗时除以本文壁钟耗时。

主表保留原五个小/中算例，增加三个 ffc/fjh 较大算例。lil 单列附录的额外枚举测试；不将其加入已构建/筛选模型集合。

| Net | C3/C4 | Assignments | Orbits (Γ) | This work / s | CASM (Γ) / s | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| ffc | 2/2 | 128 | 42 | 0.004502 | 0.0002874 | 0.06× |
| ffc | 3/3 | 2,187 | 513 | 0.004868 | 0.003095 | 0.64× |
| ffc | 5/5 | 78,125 | 14,250 | 0.01482 | 0.09148 | 6.17× |
| ffc | 5/15 | 2,109,375 | 381,375 | 0.505 | 2.753 | 5.45× |
| fjh | 1/2 | 4,096 | 394 | 0.04473 | 0.006021 | 0.13× |
| fjh | 2/1 | 65,536 | 4,756 | 0.04911 | 0.09181 | 1.87× |
| fjh | 1/3 | 531,441 | 36,369 | 0.1057 | 0.6544 | 6.19× |
| fjh | 1/4 | 16,777,216 | 1,078,816 | 2.769 | 20.94 | 7.56× |


**性能结果正文：**

The method was faster than the matched CASM workflow in five of the eight ffc and fjh cases in Table tab:casm. For fjh with |C3|=1 and |C4|=4, the 16,777,216 assignments yielded 1,078,816 orbits. Enumeration took 2.769 s, compared with 20.935 s for CASM, giving a 7.56-fold speedup. CASM was faster for the three smallest assignment spaces. Relative performance also varied within a topology: the speedup for ffc changed from 6.17 to 5.45 as the pools increased from 5/5 to 5/15.

**译文：** 在表 tab:casm 的八个 ffc 和 fjh 算例中，本文方法在五例上快于采用匹配条件的 CASM 流程。对于 |C3|=1、|C4|=4 的 fjh，16,777,216个指派产生1,078,816个轨道。本文枚举耗时2.769秒，CASM为20.935秒，对应7.56倍加速。CASM在原始指派空间最小的三个算例上更快。相对性能在同一拓扑内也会变化：ffc 的标签池由5/5增至5/15时，加速比从6.17变为5.45。

For the largest fjh case, peak process memory was 887.5 MiB for this work and 720.1 MiB for CASM. The completed-prefix cache had no hits in these benchmarks. Timing ranges, native-core times, memory measurements, and a CASM count-only comparison under Γ are reported in Supplementary Section app:enum_benchmark.

**译文：** 在最大的 fjh 算例中，本文与 CASM 的进程峰值内存分别为887.5 MiB和720.1 MiB。已完成前缀缓存在这些基准中没有命中。计时范围、原生核时间、内存测量以及 Γ 下 CASM 仅计数的对照列于补充材料 app:enum_benchmark。

## 11. 新增基准附录与原生 F 计数说明

**协议正文：**

Both workflows used Python 3.11.1 with OMP, OpenBLAS, and MKL thread counts set to one. CPU affinity was not pinned. The implementation tested here used the C++ DFS backend and the path-compressed representative Trie, with completed-prefix caching enabled. The CASM comparison called its occupation counter and canonical-configuration test with the matched subgroup, and constructed the same output fields. Count-only timings under Γ were measured separately. All exhaustive attempts finished within the external resource limits; their internal output budgets and task timeouts were disabled.

**译文：** 两种流程均使用 Python 3.11.1，并将 OMP、OpenBLAS 和 MKL 线程数设为1，未固定CPU亲和性。本文被测实现采用 C++ DFS 后端和路径压缩代表 Trie，并开启已完成前缀缓存。CASM 对照调用占位计数器及指定匹配子群的规范构型测试，构建相同的输出字段。Γ 下的仅计数时间另行测量。所有穷尽测试均在外部资源限制内完成，内部输出预算和任务超时均已关闭。

**补充原生 F 的一段：**

CASM's native full factor group F also includes orientation-reversing operations. Under F, the ffc orbit counts matched those under Γ for the tested pools, whereas the fjh counts were smaller. For fjh 1/4, the counts were 577,456 under F and 1,078,816 under Γ. These counts describe the effect of the selected equivalence group; the performance comparison in Table tab:casm uses Γ on both sides.

**译文：** CASM 的原生完整因子群 F 还包含反向操作。在被测标签池下，ffc 在 F 与 Γ 下的轨道计数相同，而 fjh 在 F 下的计数较少。对于 fjh 1/4，F 下为577,456个轨道，Γ 下为1,078,816个。这些计数说明所选等价群的影响；表 tab:casm 的性能对照两侧均使用 Γ。

**拟附表：** 原生 F 与 Γ 的计数列见下表；九例的时间 min–max、RSS、核耗时、CASM 同群仅计数与五例旧版/新版时间，直接使用已核对报告中的各表。DFS状态和CASM完整候选是不同阶段的计数，不计算它们的“总剪枝率”。

| Net | C3/C4 | Orbits (Γ) | CASM orbits (F) |
|---|---:|---:|---:|
| ffc | 2/2 | 42 | 42 |
| ffc | 3/3 | 513 | 513 |
| ffc | 5/5 | 14,250 | 14,250 |
| ffc | 5/15 | 381,375 | 381,375 |
| fjh | 1/2 | 394 | 265 |
| fjh | 2/1 | 4,756 | 2,476 |
| fjh | 1/3 | 36,369 | 20,655 |
| fjh | 1/4 | 1,078,816 | 577,456 |
| lil | 4/3 | 1,328,832 | 667,872 |


**lil 补充说明：**

An additional lil enumeration case used |C3|=4 and |C4|=3 with |Γ|=4. Its 5,308,416 assignments yielded 1,328,832 orbits. The wall times were 1.952 s for this work and 6.817 s for the matched CASM workflow. This case was used for enumeration only.

**译文：** 额外的 lil 枚举算例采用 |C3|=4、|C4|=3，且 |Γ|=4。5,308,416个指派产生1,328,832个轨道。本文和匹配 CASM 流程的壁钟耗时分别为1.952秒和6.817秒。该算例仅用于枚举。

lil 的原始位点群作用随本轮源码/输入快照保存；本轮不凭空补写其未经独立审计的空间群符号或生成元表。

## 12. vba 正文补充，保留原表配额数据

At the matched output sizes of 41,984 and 60,549, the C++ ordinary schedule returned the same ordered representatives and stabilizer metadata as the original Python schedule. Median wall times decreased from 2.248 to 0.256 s for two labels and from 3.367 to 0.281 s for three labels. Both ordinary outputs retained the stratum counts in Table tab:vba_bounded and omitted the trivial-stabilizer stratum. The quota results and task times in that table were obtained with the separate Python quota schedule.

**译文：** 在匹配输出数量41,984和60,549下，C++普通调度与原Python调度返回了相同的有序代表及稳定子元数据。两类标签的壁钟耗时中位数从2.248降为0.256秒，三类标签从3.367降为0.281秒。两组普通输出均保留了表 tab:vba_bounded 中的层计数，且均未包含平凡稳定子层。该表的配额结果和任务耗时由独立的Python配额调度得到。

**附录补充限制：** This before-and-after comparison used output caps with task timeouts disabled.

**译文：** 该优化前后对照采用输出数量上限，并关闭任务超时。

**理由：** 分开回答相同数量下的速度与各层候选分布。原 vba 表计数、括号内组成向量数、配额参数和层时间不改，不将本轮未重测的组成向量数写成已重新核验。原caption可维持，仅将 Layer time 统一为 Quota task time（配额任务耗时），如需避免改表布局则在正文完成说明即可。

## 13. 附录模板段、caption与模板池清单

**附录首段替换为：**

Table tab:building_units lists the ten node templates used for model construction and screening, together with the fixed linker template. The adsorption set contains all ten node types. Enumeration timing benchmarks additionally used larger degree-compatible label pools, whose membership is listed in Supplementary Section app:enum_benchmark.

**译文：** 表 tab:building_units 列出模型构建和筛选使用的十种节点模板，以及固定连接体模板。吸附筛选集合涉及全部十种节点类型。枚举计时基准另外采用了更大的度数相容标签池，具体成员列于补充材料 app:enum_benchmark。

**caption替换为：** Node and linker templates used for model construction and screening.

**译文：** 模型构建和筛选使用的节点与连接体模板。

**计时池说明：**

For each degree, timing pools were selected as prefixes of the available template filenames in lexicographic order. Pool cardinalities are specified in Table tab:casm; the ordered master lists are given below. The same pools were supplied to both workflows. These benchmarks enumerate assignments without constructing atomic models.

**译文：** 对每种度数，计时标签池均按可用模板文件名的字典序取前若干项。池大小列于表 tab:casm，下列主列表给出顺序。两种流程使用相同标签池。这些基准枚举指派，不构建原子模型。

**C3 顺序主列表：**

`3c_B_Ch`, `3c_Cu_imadazoline_Ch`, `3c_pyrrole_Ch`, `3c_triazine_Ch`, `3c_uranyl`。

**C4 顺序主列表：**

`4c_1anC_Ch`, `4c_4B_fused_Ch`, `4c_B_Ch`, `4c_Cd_1_Ch`, `4c_Cu_1_Ch`, `4c_In_1_Ch`, `4c_Ni_1_Ch`, `4c_V_1_Ch`, `4c_Zn_tet`, `4c_bicyclooctane_Ch`, `4c_ethene_Ch`, `4c_porphyrin_Ch`, `4c_pyrene_Ch`, `4c_pyrene_NU_Ch`, `4c_tetrathiafulvalene_Ch`。



原模板表各行 Role 仍可保留已证实的用途，不按所有新池一并扩展结构验证结论。

## 14. 结论第一段

In exhaustive runs, subgroup fixed-point searches with canonical prefix comparisons and coverage reuse return one representative per Γ-orbit. Under matched groups and output requirements, the tested implementation reached a 7.56-fold speedup over the CASM workflow for the largest fjh case, while CASM was faster for the smallest assignment spaces. Selected representatives were used to construct pre-relaxation periodic models for subsequent screening.

**译文：** 在穷尽运行中，结合规范前缀比较与覆盖复用的子群固定点搜索为每个 Γ 轨道返回一个代表。在相同作用群和输出要求下，被测实现在最大的 fjh 算例上相对 CASM 流程达到7.56倍加速，而CASM在最小的指派空间上更快。选定代表用于构建弛豫前周期模型，供后续筛选使用。

**理由：** 结论与新主表同群性能主题一致。原 Γ/F 群选择的计数讨论已移附录；结论第二段的几何审计/筛选及弛豫边界保留。

## 15. 数据可用性与最终同步检查

新增结果当前只核验了本地文件，尚未核验其已出现在公开仓库。因此本轮不把这些新增材料写成已公开。

完成数据包同步并核验公开版本后，现有 Data Availability 第二句可拟替换为：

The package includes representative and audited CIF models, enumeration benchmark scripts and raw run records, exact node-label pools and software-version manifests, source tables for the CASM comparison and DFT energies, CO2/N2 screening data, bounded-vba records, and the X-site fusion audit.

**译文：** 数据包包括代表性及审计用CIF模型、枚举基准脚本和原始运行记录、具体节点标签池及软件版本清单、CASM对照与DFT能量源表、CO2/N2筛选数据、受限vba记录，以及X位点融合审计。

英文与中文辅助注释逐处同步；实际落稿后核对：

- `alg:enum` 是普通优化算法，`alg:quota` 是配额算法，`app:quota_proof` 证明范围明确。
- 两类 Trie、规范前缀比较与已完成区域覆盖始终分开。
- 原始空间 |X|、轨道数 |X/Γ|、稳定子阶数和ρ对应现有定义。
- 主表8例、附表lil1例；主表“5/8更快”与全部九例“6/9更快”不混写。
- 本轮新标签池与原72模型审计池明确分开；构建与筛选结果不从新计时池外推。
- 计时为enumeration wall time，RSS为进程峰值；不声称全流程CIF构建提速或多项式复杂度。
- 原生F计数讨论、所有交叉引用与表格caption同步。
- 主稿英文和中文同步后，完整编译并检查算法/性能表版宽；再重新生成无注释match.tex与交付包。
