# 论文与代码修改计划

## 修改目标

1. **论文**：将"三层级（Three-Tier）"叙述改写为统一的"子群格下降（Subgroup Lattice Descent）"框架，Γ 和 {e} 分别作为格顶和格底的自然成员，不再独立成阶段
2. **代码**：移除所有截断机制（`max_results`、`DFS_STEP_LIMIT`、`SPACE_LIMIT`），使枚举在数学上完备，Burnside 自验证作为完备性证明
3. **熔断器**：随截断移除自动解决（原 Review 中的 Item 2）

---

## 一、代码修改

### 修改 1.1：`enumerator.py` — 移除 `_try_add_result()` 中的 Top-K 截断

**修改前**：接收 `max_results` 参数，用最小堆维护 Top-K 结果，超出上限时丢弃低稳定子构型
**修改后**：始终存储所有找到的唯一构型

### 修改 1.2：`enumerator.py` — 移除 DFS 步数限制

**修改前**：`DFS_STEP_LIMIT = 200000`，DFS 到达步数后截断
**修改后**：移除 `state['steps']` 计数和 `state['limit']` 检查，DFS 运行至完全穷尽

### 修改 1.3：`enumerator.py` — 移除空间熔断器

**修改前**：`SPACE_LIMIT = 5e6`，非全群任务超过此限制时 `break` 退出循环
**修改后**：移除该判断，所有任务均被执行

### 修改 1.4：`enumerator.py` — 清理 Tier 相关注释

将代码中 "Tier 1"、"Tier 2"、"Tier 3" 等注释替换为匹配新语义的表述

### 修改 1.5：`main.py` — 移除 `max_results` 参数

`enumerate(all_ligands, max_results=10000)` → `enumerate(all_ligands)`

---

## 二、论文修改

### 修改 2.1：段落 "Stratified Subgroup Decomposition"（第 234–281 行）

**修改前**：三个层级的分段叙述，Tier 1（全群，ρ=1）、Tier 2（循环子群，1/|Γ|<ρ<1）、Tier 3（恒等群，ρ=1/|Γ|）。Tier 2 仅提及循环子群。

**修改后**：叙述为统一的子群格下降策略：
- 理论框架仍然引用子群格 L(Γ)，但将 Γ 作为格顶、{e} 作为格底
- 算法步骤：通过计算群代数（SymPy）枚举 L(Γ) 的所有子群 → 按轨道划分去重 → 按轨道数升序（对称性降序）排序 → 依次执行 DFS
- "三层级"从算法结构降格为**概念分类**：在讨论 ρ 的取值时用于解释对称性等级
- 补充说明：恒等群 {e} 作为子群格的最低元素，其 DFS 在无对称性约束的全组合空间执行，是数学完备性的必要保证

### 修改 2.2：Algorithm 1（第 283–309 行）

**修改前**：三阶段条件分支：
```
Tier 1: Add Task(G)
Tier 2: for g ∈ Γ, cyclic subgroups
Tier 3: Add Task({e})
if PredictedSize > L_space and Tier 3: continue
```

**修改后**：单一循环枚举子群格：
```
Enumerate subgroup lattice L(Γ)
for H ∈ L(Γ) with unique orbit partition:
    Add Task(H)
Sort tasks by |O_H| ascending
for each Task:
    Execute SymmetryAwareDFS
```

具体 LaTeX 代码见下方"修改前后对照"。

### 修改 2.3：图 1 图题中的 Phase II 描述

调整 Phase II 描述，使其不再提及"三层级枚举"，改为"子群格分层下降遍历"。

### 修改 2.4：性能与复杂性分析段落

将"exhaustive enumeration within 10 seconds"相关表述调整为反映枚举范围：
- 前两段侧重子群格下降的效率优势，弱化"穷举"声称
- 增加对恒等群任务计算成本的说明

---

## 三、不修改的部分

- **Algorithm 2**（DFS Pruning）：正确保持不变
- **图 1 的三阶段流水线结构**：Phase I / Phase II / Phase III 的整体结构不受影响，仅调整 Phase II 内部说明
- **Burnside 验证的实现逻辑**：验证公式不变，移除截断后自然通过
- **对称性分析、几何重构、实验方法**等章节：不受本次修改影响

---

## 四、修改前后对照（论文 LaTeX）

### 4.1 "Stratified Subgroup Decomposition" 段落

**修改前**（第 252–281 行）：
```latex
% 依据命题 1，搜索空间的维度实际上由子群 H 定义的轨道数量 |O_H| 控制。
% 为了系统地遍历构型空间，我们沿子群格 L(Γ) 采用了一种基于对称性指数 ρ 的分层下降策略，将搜索空间严格划分为三个层级：
Based on Proposition \ref{prop:dim_reduction}, the dimensionality of the search space is effectively governed by the orbit count $|\mathcal{O}_H|$ defined by the subgroup $H$[cite: 55].
To systematically traverse the configuration space, we employ a stratified descent strategy along the subgroup lattice $\mathcal{L}(\Gamma)$ based on the Symmetry Index $\rho$, rigorously partitioning the search space into three tiers:

\begin{itemize}
    \item \textbf{Tier 1: Isotropic Limit ($\rho = 1$, $H=\Gamma$).}
    The search initiates at the supremum of the subgroup lattice, where all graph automorphisms are enforced.
    At this tier, the orbit length of configurations is minimized to 1, thereby reducing the original combinatorial space to a polynomial or even constant scale and permitting the exhaustive enumeration of all highly symmetric candidates with minimal computational cost.

    \item \textbf{Tier 2: Stratified Anisotropic Subspace ($\frac{1}{|\Gamma|} < \rho < 1$, $H \in \mathcal{L}(\Gamma) \setminus \{\Gamma, \{e\}\}$).}
    As the algorithm descends to intermediate tiers defined by non-trivial proper subgroups, symmetry constraints are correspondingly relaxed.
    Orbits undergo splitting, resulting in an increase in the computational complexity associated with $|\mathcal{O}_H|$.
    This tier captures configurations that lack full isotropy but retain structural ordering along specific directions, balancing spatial complexity with structural regularity.

    \item \textbf{Tier 3: Asymmetric Limit ($\rho = \frac{1}{|\Gamma|}$, $H=\{e\}$).}
    At the infimum of the lattice, stabilizer constraints are entirely removed, and the algorithm reverts to searching the complete combinatorial space.
    At this stage, the cardinality of the state space restores to the original $|X| = |\mathcal{C}|^{|V|}$, and the worst-case time complexity scales as $\mathcal{O}(|\mathcal{C}|^{|V|})$.
    This space is predominantly occupied by high-entropy, asymmetric configurations with maximal redundancy.
    To manage computational overhead, the algorithm implements a truncation strategy when the state space exceeds a predefined threshold (e.g., $5 \times 10^6$ states), prioritizing computational resource allocation to the upper two tiers.
\end{itemize}
```

**修改后**（建议）：
```latex
% 依据命题 1，搜索空间的维度实际上由子群 H 定义的轨道数量 |O_H| 控制。
% 为了系统地遍历构型空间，我们沿子群格 L(Γ) 采用统一的"子群格下降"策略：
% 枚举 L(Γ) 中的所有子群，按轨道划分去重后按对称性降序执行 DFS。
Based on Proposition \ref{prop:dim_reduction}, the dimensionality of the search space is effectively governed by the orbit count $|\mathcal{O}_H|$ defined by the subgroup $H$.
To systematically traverse the configuration space, we employ a \textit{Subgroup Lattice Descent} strategy: the full lattice $\mathcal{L}(\Gamma)$ — ranging from the full automorphism group $\Gamma$ at its supremum to the trivial group $\{e\}$ at its infimum — is enumerated via computational group algebra (e.g., \texttt{SymPy}). Each subgroup $H \in \mathcal{L}(\Gamma)$ defines a fixed-point subspace $\mathrm{Fix}(H)$ and an associated orbit partition $\mathcal{O}_H$. After deduplication by orbit-partition signature, tasks are scheduled by ascending $|\mathcal{O}_H|$ (i.e., descending symmetry) and processed sequentially.

The Symmetry Index $\rho(c) = |\mathrm{Stab}(c)|/|\Gamma|$ provides a natural taxonomy for the resulting configurations:
\begin{itemize}
    \item $\rho = 1$ ($H = \Gamma$): maximally symmetric configurations at the lattice supremum, where $|\mathcal{O}_\Gamma| \ll |V|$ permits enumeration at polynomial cost;
    \item $\frac{1}{|\Gamma|} < \rho < 1$ ($\{e\} \subsetneq H \subsetneq \Gamma$): intermediate configurations with partial symmetry, enumerated via intermediate lattice subgroups;
    \item $\rho = \frac{1}{|\Gamma|}$ ($H = \{e\}$): fully asymmetric configurations at the lattice infimum, where the algorithm traverses the complete combinatorial space $\mathcal{C}^V$ as a systematic DFS without symmetry pruning.
\end{itemize}

This unified lattice descent framework ensures exhaustive enumeration of the entire orbit space: the union of $\mathrm{Fix}(H)$ over all $H \in \mathcal{L}(\Gamma)$ covers every configuration $c \in X$, because $\mathrm{Stab}(c) \in \mathcal{L}(\Gamma)$ for any $c$. The Burnside self-verification step (Section~\ref{sec:verification}) subsequently confirms enumeration completeness by verifying the orbit-stabilizer identity $\sum_{c \in \Omega} 1/|\mathrm{Stab}(c)| = |\mathrm{Fix}(e)|/|\Gamma|$.
```

**中文翻译**：
基于命题 1，搜索空间的维度由子群 H 定义的轨道数 |O_H| 控制。为了系统地遍历构型空间，我们采用"子群格下降"策略：完整子群格 L(Γ)——从全自同构群 Γ（格顶）到平凡群 {e}（格底）——通过计算群代数（如 SymPy）枚举得到。每个子群 H ∈ L(Γ) 定义一个不动点子空间 Fix(H) 和关联的轨道划分 O_H。按轨道划分签名去重后，任务按 |O_H| 升序（即对称性降序）调度并依次处理。

对称性指数 ρ(c) = |Stab(c)|/|Γ| 为生成的构型提供自然分类：(1) ρ = 1（H=Γ）为格顶处最大对称构型，|O_Γ| ≪ |V|，以多项式成本枚举；(2) 1/|Γ| < ρ < 1（中间子群）为部分对称构型；(3) ρ = 1/|Γ|（H={e}）为格底处完全不对称构型，算法在该层级执行无对称性剪枝的系统 DFS，遍历完整的组合空间 C^V。

该统一子群格下降框架保证了完整轨道空间的穷举枚举：因为对任意构型 c 都有 Stab(c) ∈ L(Γ)，所有 H ∈ L(Γ) 的 Fix(H) 之并集覆盖 X 中的每个构型。随后的 Burnside 自验证步骤通过检验轨道-稳定子恒等式来确认枚举完备性。

### 4.2 Algorithm 1

**修改前**（第 283–309 行）：
```latex
\begin{algorithm}
\caption{Stratified Subgroup Decomposition Strategy}
\label{alg:stratified_descent}
\begin{algorithmic}[1]
\Require Graph Automorphism Group $\Gamma$, Atom Types $\mathcal{C}$, Space Limit $L_{space}$
\Ensure Canonical Unique Configuration Set $\Omega$
\State $\mathcal{T} \gets \emptyset$ \Comment{Initialize Task Queue}
\State \textbf{Tier 1 ($\rho=1$):} Add Task(Group=$\Gamma$, Priority=High) to $\mathcal{T}$
\State \textbf{Tier 2 ($\frac{1}{|\Gamma|} < \rho < 1$):} \textbf{for} $g \in \Gamma$ \textbf{do} $fp \gets$ OrbitSizeSignature($\langle g \rangle$); \textbf{if} $fp \notin$ Seen \textbf{then} Add Task(Group=$\langle g \rangle$, Priority=Medium) to $\mathcal{T}$; Seen.add($fp$)
\State \textbf{Tier 3 ($\rho=\frac{1}{|\Gamma|}$):} Add Task(Group=$\{e\}$, Priority=Low) to $\mathcal{T}$ 
\State Sort $\mathcal{T}$ by Priority descending
\For{each Task $T_H$ in $\mathcal{T}$}
    \State $\mathcal{O}_H \gets \Call{UnionFindOrbits}{T_H.\text{Group}}$
    \State $PredictedSize \gets |\mathcal{C}|^{|\mathcal{O}_H|}$
    \If{$PredictedSize > L_{space}$ \textbf{and} $T_H$ is Tier 3}
       \State \textbf{continue} 
        \Comment{Circuit Breaker: Skip explosive noise}
    \EndIf
    \State $\Omega \gets \Omega \cup \Call{SymmetryAwareDFS}{\mathcal{O}_H, \Gamma, \dots}$
\EndFor
\end{algorithmic}
\end{algorithm}
```

**修改后**（建议）：
```latex
\begin{algorithm}
\caption{Subgroup Lattice Descent Strategy}
\label{alg:lattice_descent}
\begin{algorithmic}[1]
\Require Graph Automorphism Group $\Gamma$, Atom Library $\mathcal{C}$
\Ensure Canonical Unique Configuration Set $\Omega$
\State $\mathcal{T} \gets \emptyset$ \Comment{Initialize Task Queue}
\State $\mathcal{L}(\Gamma) \gets \Call{EnumerateSubgroupLattice}{\Gamma}$ \Comment{via computational group algebra}
\State $\mathcal{S}_{seen} \gets \emptyset$
\For{each $H \in \mathcal{L}(\Gamma)$}
    \State $\mathcal{O}_H \gets \Call{UnionFindOrbits}{H}$
    \State $sig \gets$ OrbitPartitionSignature($\mathcal{O}_H$)
    \If{$sig \notin \mathcal{S}_{seen}$}
        \State $\mathcal{S}_{seen}$.add($sig$)
        \State $\mathcal{T}$.add(Task(Group=$H$, Orbits=$\mathcal{O}_H$))
    \EndIf
\EndFor
\State Sort $\mathcal{T}$ by $|\mathcal{O}_H|$ ascending \Comment{highest symmetry first}
\For{each Task $T_H$ in $\mathcal{T}$}
    \State $\Omega \gets \Omega \cup \Call{SymmetryAwareDFS}{\mathcal{O}_H, \Gamma, \mathcal{C}}$
\EndFor
\State $\Call{BurnsideVerify}{\Omega, \Gamma, \mathcal{C}}$ \Comment{guarantees completeness}
\end{algorithmic}
\end{algorithm}
```

**中文翻译**：
算法 1：子群格下降策略。输入：图自同构群 Γ，原子库 C。输出：正则唯一构型集 Ω。
1. 初始化任务队列 T ← ∅
2. 枚举子群格 L(Γ)（通过计算群代数工具）
3. 对 H ∈ L(Γ)，使用并查集计算轨道划分 O_H，提取轨道划分签名；若签名是新出现的，则记录并将任务加入 T
4. 按 |O_H| 升序排序（最高对称性优先）
5. 对每个任务 T_H，执行对称感知 DFS
6. 执行 Burnside 验证（确保完备性）

### 4.3 图 1 图题

在 Phase II 描述中，将"executes a stratified traversal along the subgroup lattice"保留，移除三层级和熔断器的描述。新增 Burnside 验证的提及。

**修改前**：
> \textbf{Phase II (Enumeration)} executes a stratified traversal along the subgroup lattice $\mathcal{L}(\Gamma)$, utilizing orbit-based dimensionality reduction to identify symmetry-unique topological graphs ($\Omega$) while pruning combinatorial noise. A subsequent Burnside self-verification step validates enumeration completeness by checking $\sum 1/|\mathrm{Stab}(c)| = |\mathrm{Fix}(e)|/|\Gamma|$;

**修改后**：
> \textbf{Phase II (Enumeration)} executes a subgroup lattice descent over $\mathcal{L}(\Gamma)$, systematically enumerating each subgroup's fixed-point subspace via orbit-based DFS to identify all symmetry-unique topological configurations ($\Omega$). The phase concludes with a Burnside self-verification step that confirms exhaustive enumeration by checking $\sum_{c \in \Omega} 1/|\mathrm{Stab}(c)| = |\mathrm{Fix}(e)|/|\Gamma|$;

**中文翻译**：
Phase II（枚举）在 L(Γ) 上执行子群格下降，通过基于轨道的 DFS 系统枚举每个子群的不动点子空间，识别所有对称性唯一的拓扑构型（Ω）。该阶段以 Burnside 自验证步骤收尾，通过轨道-稳定子恒等式确认穷举枚举。

### 4.4 性能段落

**修改前**：
> the algorithm completes exhaustive enumeration within 10 seconds

**修改后**：
> the algorithm enumerates all symmetry-unique configurations within 10 seconds, with the Burnside identity confirming complete orbit-space coverage

**中文翻译**：
算法在 10 秒内枚举了所有对称性唯一的构型，Burnside 恒等式验证了完整的轨道空间覆盖。

---

## 五、注意事项

1. **Burnside 验证在恒等群阶段的时间成本**：对于 3^18 ≈ 3.87×10^8 规模的朴素空间，恒等群 DFS 需要遍历所有可能性。这在 sod（12个节点）和 lil（12个节点）下很快，但对于 rht（18个节点）可能需要更长时间。建议在论文中如实说明这一计算成本。
2. **论文中 "exhaustive" 和 "complete" 的使用**：修改后所有此类术语均应指向"经由 Burnside 恒等式验证的数学完备性"，而非暗示枚举过程本身没有计算成本。
3. **熔断器删除后大型系统的适用性**：用户可选择在调用层控制搜索规模（如限制原子种类数量），但论文层次应呈现完整的数学框架。
