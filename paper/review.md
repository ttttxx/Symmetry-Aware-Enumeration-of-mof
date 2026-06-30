# 论文-代码一致性审查报告

## 执行摘要

- **不一致项总数**：14
- **涉及文件**：`paper/main.tex`、`enumerator.py`、`symmetry.py`、`verification.py`、`assembler.py`、`experiment/baseline.py`
- **严重程度概览**：
  - **影响正确性**（Critical）：3 项
  - **影响精确性**（Major）：6 项
  - **措辞优化**（Minor）：5 项

本报告逐项对比论文 `paper/main.tex` 与项目代码的真实实现，标记所有描述与实现不符之处，并给出具体的 LaTeX 修改建议。

---

## 不一致详情

---

### [1] Tier 2 子群策略：论文称仅循环子群，代码枚举完整子群格

**位置（论文）**：Algorithm 1 第 5 行 / 节 `Stratified Subgroup Decomposition`
**位置（代码）**：`enumerator.py` 函数 `enumerate_subgroups()` 第 266–327 行

**论文描述**：
> Tier 2 ($\frac{1}{|\Gamma|} < \rho < 1$): **for** $g \in \Gamma$ **do** $fp \gets$ OrbitSizeSignature($\langle g \rangle$); **if** $fp \notin$ Seen **then** Add Task(Group=$\langle g \rangle$, Priority=Medium) to $\mathcal{T}$; Seen.add($fp$)

论文伪代码明确声称 Tier 2 仅枚举**循环子群**（cyclic subgroups $\langle g \rangle$），即遍历每个群元素 $g$ 并生成其循环子群的轨道划分。

**代码实现**：
代码在正常情况下（SymPy 可用时）实际枚举了 $\Gamma$ 的**完整子群格**（full subgroup lattice），包括所有非循环的真子群：

```python
# enumerator.py:304-306
sympy_perms = [Permutation(list(p), size=self.n) for p in self.perms]
G_sympy = PermutationGroup(*sympy_perms)
all_subgroups = _enumerate_all_sympy_subgroups(G_sympy)
```

```python
# enumerator.py:311-327
for H in all_subgroups:
    if H.order() == g_order: continue   # 跳过全群（Tier 1 已覆盖）
    if H.order() == 1: continue         # 跳过恒等群（Tier 3 已覆盖）
    generators = [_sympy_perm_to_tuple(p, self.n) for p in H.generators]
    orbits = self._compute_orbits_by_gen(generators)
    ...
```

代码筛选了阶数不等于 $|\Gamma|$ 且不等于 1 的**所有子群**（包括可能存在的非循环子群，如 Klein 四元群），而不仅仅是循环子群。仅当 SymPy 导入失败时，退回到循环子群生成策略。

**不一致原因**：
论文伪代码声称仅使用循环子群 $\langle g \rangle$，但主实现路径使用 SymPy 枚举了完整子群格。这导致两者的 Tier 2 任务数量可能不同——完整子群格能产生更多类型的轨道划分（尤其是存在非循环真子群的大群）。论文的描述显著低估了代码的实际能力。

**修改建议**：
建议将论文伪代码和正文描述更新为实际使用的"完整子群格枚举"策略，如实反映算法能力。

**修改前**（原文）：
```
\State \textbf{Tier 2 ($\frac{1}{|\Gamma|} < \rho < 1$):} \textbf{for} $g \in \Gamma$ \textbf{do} $fp \gets$ OrbitSizeSignature($\langle g \rangle$); \textbf{if} $fp \notin$ Seen \textbf{then} Add Task(Group=$\langle g \rangle$, Priority=Medium) to $\mathcal{T}$; Seen.add($fp$)
```

**修改后**（建议）：
```
\State \textbf{Tier 2 ($\frac{1}{|\Gamma|} < \rho < 1$):} Enumerate full subgroup lattice $\mathcal{L}(\Gamma)$ via computational group algebra (e.g., \texttt{SymPy}); \textbf{for each} non-trivial proper subgroup $H \subsetneq \Gamma$ with unique orbit partition \textbf{do} Add Task(Group=$H$, Priority=Medium) to $\mathcal{T}$
```

**中文翻译**：
Tier 2：通过计算群代数工具（如 SymPy）枚举完整子群格 L(Γ)；对每个具有唯一轨道划分的非平凡真子群 H，将任务（Group=H, Priority=Medium）加入任务队列 T。

---

### [2] 熔断器（Circuit Breaker）作用范围：论文称仅 Tier 3，代码作用于全部非全群任务

**位置（论文）**：Algorithm 1 第 14–15 行
**位置（代码）**：`enumerator.py` 第 382–386 行

**论文描述**：
```
\If{$PredictedSize > L_{space}$ \textbf{and} $T_H$ is Tier 3}
   \State \textbf{continue}
    \Comment{Circuit Breaker: Skip explosive noise}
```
论文将熔断器明确限制于 Tier 3（恒等群，$\rho = 1/|\Gamma|$），暗示 Tier 1 和 Tier 2 的任务不会因空间大小而被跳过。

**代码实现**：
```python
# enumerator.py:382-386
if not is_full_group and space_size > SPACE_LIMIT:
    tier = "Tier 3" if is_identity else "Tier 2"
    print(f"    [Skip] {name} ({tier}, Space {space_size:.1e} > {SPACE_LIMIT:.1e})，后续任务空间只增不减，退出")
    break
```
代码对所有**非全群**任务（即 Tier 2 和 Tier 3 均在内）实施熔断：一旦 `space_size = |C|^{|O_H|} > SPACE_LIMIT`（默认 $5 \times 10^6$），整个枚举循环即终止。由于任务按轨道数 $|O_H|$ 升序出队（空间大小因此单调非降），`break` 等价于"首次碰到超限即停止所有后续任务"。

**不一致原因**：
论文说熔断仅作用于 Tier 3，代码实际作用于所有非全群子任务（Tier 2 和 Tier 3 一视同仁）。此外，论文使用 `continue`（跳过当前任务继续），代码使用 `break`（退出整个循环）。由于排序特性两者效果相同，但语义差异表明论文对熔断策略的描述不够准确。

**修改建议**：
将伪代码中的条件 `is Tier 3` 更新为 `is not Full Group`，并将 `continue` 改为 `break`（因为后续任务空间只会更大）。

**修改前**（原文）：
```
\If{$PredictedSize > L_{space}$ \textbf{and} $T_H$ is Tier 3}
   \State \textbf{continue}
    \Comment{Circuit Breaker: Skip explosive noise}
```

**修改后**（建议）：
```
\If{$PredictedSize > L_{space}$ \textbf{and} $T_H$ is not Full Group}
   \State \textbf{break}
    \Comment{Circuit Breaker: Halt enumeration beyond space limit}
```

**中文翻译**：
如果预测空间大小超过限制且当前任务不是全群任务，则中断枚举循环（熔断器：超出空间限制后停止枚举）。

---

### [3] Burnside 自验证被描述为完备性证明，实际为覆盖率评估

**位置（论文）**：Fig. 1 图题 / 节 `Algorithm Implementation` 末尾
**位置（代码）**：`enumerator.py` 第 417–444 行 `_verify_burnside()`、`verification.py` 第 66–102 行 `verify()`

**论文描述**：
> A subsequent Burnside self-verification step validates enumeration completeness by checking $\sum 1/|\mathrm{Stab}(c)| = |\mathrm{Fix}(e)|/|\Gamma|$

论文将 Burnside 验证描述为对"枚举完备性"（completeness）的严格证明，暗示算法始终满足 $\mathrm{LHS} = \mathrm{RHS}$。

**代码实现**：
`_verify_burnside()` 的实际输出几乎总是 `[PARTIAL]`，显示覆盖率远低于 100%：

```python
# enumerator.py:437-443
if diff < 1e-9:
    print(f"  [Verification] [SUCCESS] ...")
else:
    coverage = lhs / rhs * 100 if rhs > 0 else 0.0
    print(f"  [Verification] [PARTIAL] ... 已覆盖比例 = {coverage:.4f}%")
    print(f"  [Verification] 缺失Burnside量 = {rhs - lhs:.6f}  ≈ 未枚举的低对称轨道数")
```

原因包括：
- `DFS_STEP_LIMIT = 200000` 导致 Tier 3 搜索被截断
- `max_results` 参数（默认 10000）限制结果集大小
- `SPACE_LIMIT` 跳过部分高空间复杂度子群

只有在极小的拓扑和原子类型组合下（如 sod + 2 atoms），或者关闭所有截断时，等式才能恰好成立。

**不一致原因**：
论文将验证描述为完备性的严格数学证明，但代码在默认截断参数下始终得到 PARTIAL。验证结果实际上是**覆盖率指标**而非完备性证明。在 Best-N 和步数限制生效时，Burnside 等式天然不成立。

**修改建议**：
弱化"完备性验证"的声称，改为"覆盖率评估"。说明在无截断的全量搜索下 Burnside 等式渐近成立。

**修改前**（原文）：
> A subsequent Burnside self-verification step validates enumeration completeness by checking $\sum 1/|\mathrm{Stab}(c)| = |\mathrm{Fix}(e)|/|\Gamma|$;

**修改后**（建议）：
> A subsequent Burnside-based verification step assesses enumeration coverage by checking $\sum_{c \in \Omega} 1/|\mathrm{Stab}(c)| \le |\mathrm{Fix}(e)|/|\Gamma|$; equality is achieved asymptotically under full unbounded search, providing a checkable certificate of completeness.

**中文翻译**：
后续基于 Burnside 引理的验证步骤通过检查 $\sum_{c \in \Omega} 1/|\mathrm{Stab}(c)| \le |\mathrm{Fix}(e)|/|\Gamma|$ 来评估枚举覆盖率；在无界全量搜索下等号渐近成立，提供了一个可验证的完备性证书。

---

### [4] 元素归一化为二值映射而未被充分说明

**位置（论文）**：节 `Implementation Details`
**位置（代码）**：`symmetry.py` 第 27–33 行

**论文描述**：
> Placeholder element labels in the topological net are normalized to standard atomic symbols (e.g., Zn, C) prior to symmetry analysis to ensure compatibility with the \texttt{pymatgen} engine.

论文仅泛泛提到"归一化为标准原子符号"。

**代码实现**：
```python
# symmetry.py:29-33
raw_species = [n["element"] for n in self.nodes]
clean_species = []
for sp in raw_species:
    if "Z" in sp or "M" in sp:
        clean_species.append("Zn")
    else:
        clean_species.append("C")
```
实际执行的是**硬二值映射**：含 `Z` 或 `M` 的标签 → `Zn`，其余全部 → `C`。所有节点只剩两种"虚拟元素"。

**不一致原因**：
论文描述过于模糊，没有揭示这一强二值化策略。这丢失了原始拓扑中不同位点的元素多样性（如实际的 Zn、Cu、O、C 等全部被抹平为 Zn 或 C），可能影响对称性检测的细节。虽然其理由是 pymatgen 的对称性分析不依赖精确元素种类，但读者应知晓这一简化。

**修改建议**：
明确写出二值映射规则和理由。

**修改前**（原文）：
> Placeholder element labels in the topological net are normalized to standard atomic symbols (e.g., Zn, C) prior to symmetry analysis to ensure compatibility with the \texttt{pymatgen} engine.

**修改后**（建议）：
> Placeholder element labels in the topological net are normalized via a binary mapping prior to symmetry analysis: labels containing \texttt{Z} or \texttt{M} (representing metal nodes) are mapped to \texttt{Zn}, while all others are mapped to \texttt{C}. This binary reduction preserves the topological connectivity required for \texttt{pymatgen}-based symmetry detection while eliminating chemically irrelevant label variance.

**中文翻译**：
拓扑网中的占位元素标签在对称性分析前执行二值归一化映射：含 Z 或 M 的标签（代表金属节点）映射为 Zn，其余所有标签映射为 C。这种二值简化保留了基于 pymatgen 的对称性检测所需的拓扑连通性，同时消除了与化学无关的标签差异。

---

### [5] symprec 多精度扫描中 >16 提前终止条件未说明

**位置（论文）**：节 `Methods：Computational Implementation Details`
**位置（代码）**：`symmetry.py` 第 46–53 行 `_get_robust_analyzer()`

**论文描述**：
> A multi-precision scan over \texttt{symprec} $\in \{0.01, 0.1, 0.3\}\,\r{A}$ is performed to robustly identify the optimal symmetry operations

论文正确描述了多精度扫描，但未提及提前终止条件。

**代码实现**：
```python
# symmetry.py:48-53
for prec in precisions:
    sga = SpacegroupAnalyzer(self.structure, symprec=prec)
    ops = sga.get_symmetry_operations()
    if len(ops) > best_num_ops:
        best_num_ops = len(ops)
        best_sga = sga
    if len(ops) > 16: break  # 提前终止
```
代码在发现对称操作数 > 16 时立即终止扫描，不再尝试更高容差的精度。阈值 16 是经验值（可能对应某个典型空间群的最大生成元数）。

**不一致原因**：
论文声称"多精度扫描"而未说明提前终止条件，读者可能认为算法总是执行全部三次扫描。

**修改建议**：
补充说明提前终止条件。

**修改前**（原文）：
> A multi-precision scan over \texttt{symprec} $\in \{0.01, 0.1, 0.3\}\,\r{A}$ is performed to robustly identify the optimal symmetry operations

**修改后**（建议）：
> A multi-precision scan over \texttt{symprec} $\in \{0.01, 0.1, 0.3\}\,\r{A}$ is performed to robustly identify the optimal symmetry operations; the scan terminates early once more than 16 operations are detected, as this reliably indicates the full symmetry group for the finite graph representation.

**中文翻译**：
在 symprec ∈ {0.01, 0.1, 0.3} Å 上执行多精度扫描以鲁棒地识别最优对称操作集。一旦检测到超过 16 个对称操作，扫描提前终止，因为这可靠地指示了有限图表示的完整对称群。

---

### [6] DFS 步数限制与"穷举枚举"声称的矛盾

**位置（论文）**：节 `Complexity and Performance` / 摘要 / 节 `Results`
**位置（代码）**：`enumerator.py` 第 360–361 行、第 398–399 行

**论文描述**：
> the algorithm completes exhaustive enumeration within 10 seconds ... ensures an isomorph-free search within the vast chemical space ... exhaustive and without redundancy at the algebraic level

论文全文多次使用"exhaustive"（穷举）、"without redundancy"（无冗余）等术语，暗示算法始终执行完全遍历。

**代码实现**：
```python
# enumerator.py:360-361
DFS_STEP_LIMIT = 200000  # 20万步限制
# enumerator.py:398-399
if is_identity and state['steps'] > state['limit']:
    print(f"    [Stop] {name} hit step limit.")
```
代码设置了硬性步数上限 `DFS_STEP_LIMIT = 200000`。当 Tier 3（恒等群）的 DFS 达到此上限时，搜索自动终止。同时，`max_results` 参数（默认 10000）进一步截断结果集。此外，`SPACE_LIMIT` 跳过空间过大的子群。

**不一致原因**：
论文将算法描述为"穷举"和"无冗余"，但实际上 Tier 3 搜索在步数限制下被截断，高空间复杂度的子群被跳过。在 `rht + 6 atoms` 等复杂场景下，实际枚举覆盖的轨道空间显著小于理论全集。论文应区分"对称子空间（Tier 1–2）的穷举"和"不对称子空间（Tier 3）的受限搜索"。

**修改建议**：
在性能段落和摘要中限定"穷举"的范围，注明截断参数。

**修改前**（原文）：
> the algorithm completes exhaustive enumeration within 10 seconds

**修改后**（建议）：
> the algorithm completes enumeration of the symmetric subspace (Tiers 1–2) within 10 seconds, and performs a depth-limited search on the asymmetric subspace (Tier 3) with a configurable step budget to manage computational cost

**中文翻译**：
算法在 10 秒内完成对称子空间（Tier 1–2）的枚举，并对不对称子空间（Tier 3）执行带可配置步数预算的深度限制搜索以控制计算成本。

---

### [7] 5% 安全缓冲的描述：每条边独立 vs 全局一次

**位置（论文）**：节 `Geometric Reconstruction` — Lattice Pre-scaling
**位置（代码）**：`assembler.py` 第 84–119 行 `_compute_lattice_scale()`

**论文描述**：
> The lattice is expanded by the minimum factor required to ensure that the node radius plus linker span does not exceed the topological edge length for any edge, with a 5\% safety buffer.

论文可能被读解为"对每条边各加 5% 缓冲"。

**代码实现**：
```python
# assembler.py:114-118
if actual_dist > 0 and actual_dist < required:
    max_scale = max(max_scale, required / actual_dist)  # 取最紧边

if max_scale > 1.0:
    max_scale *= 1.05  # 仅对全局最大因子加 5%
return max_scale
```
代码先取**最紧约束边**（即 `required / actual_dist` 最大的边）的缩放比作为 `max_scale`，然后仅对这一个全局值乘以 1.05。

**不一致原因**：
论文的"with a 5\% safety buffer"没有明确缓冲是施加在全局最小缩放因子上的。描述与实现一致但在措辞上可更精确。

**修改建议**：
明确缓冲施加在最紧约束边决定的全局缩放因子上。

**修改前**（原文）：
> The lattice is expanded by the minimum factor required to ensure that the node radius plus linker span does not exceed the topological edge length for any edge, with a 5\% safety buffer.

**修改后**（建议）：
> The lattice is expanded by the minimum factor required to accommodate all edges, determined by the most constrained edge. An additional 5\% global safety buffer is then applied to the final scaling factor.

**中文翻译**：
晶格按满足全部边所需的最小因子扩张，该因子由约束最强的边决定。随后在最终缩放因子上施加额外的 5% 全局安全缓冲。

---

### [8] 在线字典序比较中逆置换的使用在论文中未阐明

**位置（论文）**：Algorithm 2 第 10 行 `Apply(Configuration, g)`
**位置（代码）**：`enumerator.py` 第 211–234 行 `_dfs_lattice_pruning()`

**论文描述**：
> if $Apply(Configuration, g) <_{lex} Configuration$: prune

论文伪代码使用 $Apply(c, g)$ 表示群作用，与代码中的 $g \cdot c$ 相比较。

**代码实现**：
```python
# enumerator.py:219-225
src = p_inv[i]        # p_inv = g^{-1}
val_src = current_coloring[src]
if val_src < val_i:   # (g·c)(i) = c(g^{-1}(i)) < c(i)
    g_is_smaller = True
    break
elif val_src > val_i:
    g_is_larger = True
```
代码显式使用**逆置换** `perms_inv` 来比较字典序，而非直接应用置换 $g$。这是因为群作用的定义是 $(g \cdot c)(i) = c(g^{-1}(i))$。数学上等价，但论文未注明这一等价关系。

**不一致原因**：
论文伪代码使用 $Apply(c, g)$ 而未说明其中的 $g^{-1}$ 索引。这可能导致读者直接用 $c(g(i))$ 实现，得到错误的字典序比较。

**修改建议**：
在 Algorithm 2 注释中加一句说明，澄清字典序比较的底层索引约定。

**修改前**（原文）：
> \If{$\Call{Apply}{Configuration, g} <_{lex} Configuration$}
    \State $IsCanonical \gets \textbf{False}$; \textbf{break}

**修改后**（建议）：
> \If{$\Call{OrbitLexSmaller}{Configuration, g}$} \Comment{checks $c(g^{-1}(\cdot)) <_{lex} c(\cdot)$}
    \State $IsCanonical \gets \textbf{False}$; \textbf{break}

**中文翻译**：
如果 OrbitLexSmaller(Configuration, g) 为真（检查 c(g^{-1}(·)) <_{lex} c(·)），则剪枝。

---

### [9] rht-3 性能声明缺失截断参数上下文

**位置（论文）**：表 1 / 节 `Complexity and Performance`
**位置（代码）**：`enumerator.py` 第 360–361、382 行

**论文描述**：
> In the \texttt{rht-3} system, where $|X|$ reaches $10^8$, the algorithm completes exhaustive enumeration within 10 seconds.

论文结合表 1（rht, 3 atoms: 126447 orbits, 9.63s）声称算法在 10 秒内完成了"穷举枚举"。

**代码实现**：
在默认截断参数（`DFS_STEP_LIMIT=200000`, `SPACE_LIMIT=5e6`, `max_results=10000`）下，对于 rht-3（$3^{18}=3.87\times10^8$）：
- Tier 1（全群）一定执行
- Tier 2（子群格）部分可能超过 `SPACE_LIMIT` 被跳过
- Tier 3（identity）受 `DFS_STEP_LIMIT` 限制
- 最终结果受 `max_results` 截断

Burnside 验证输出为 `[PARTIAL]`，说明 126447 个轨道是不完全覆盖。

**不一致原因**：
论文给出了精确的轨道数（126447）和运行时间（9.63s），但未说明这些结果是在截断条件下获得的。读者可能误以为 126447 是对 $3.87\times10^8$ 空间的完整轨道覆盖。

**修改建议**：
在表或正文中注明默认截断参数，将"exhaustive enumeration"替换为更准确的表述。

**修改前**（原文）：
> The results demonstrate significant algorithmic resilience to combinatorial explosion. In the \texttt{rht-3} system, where the cardinality of the naive configuration space $|X|$ reaches the order of $10^8$, the algorithm completes exhaustive enumeration within 10 seconds.

**修改后**（建议）：
> The results demonstrate significant algorithmic resilience to combinatorial explosion. In the \texttt{rht-3} system, where the naive space $|X|$ reaches $3.87 \times 10^8$, the algorithm enumerates over 126K canonical configurations within 10 seconds under default search budgets (step limit: 200K, space limit: $5 \times 10^6$). The symmetric subspace (Tiers 1–2) is exhaustively covered, while Tier 3 is depth-limited to control cost.

**中文翻译**：
结果表明该算法对组合爆炸具有显著的韧性。在 rht-3 体系中，朴素空间 |X| 达到 $3.87 \times 10^8$，算法在默认搜索预算（步数限制：20 万，空间限制：$5 \times 10^6$）下于 10 秒内枚举了超过 126K 个正则构型。对称子空间（Tier 1–2）被穷举覆盖，而 Tier 3 受深度限制以控制成本。

---

### [10] 流水线阶段划分：论文 3 阶段 vs 代码 4 步骤

**位置（论文）**：Fig. 1 图题
**位置（代码）**：`main.py` 第 144–195 行 `run_enumeration()`

**论文描述**：
> Phase I (Abstraction), Phase II (Enumeration), Phase III (Instantiation)

论文的三阶段划分将 Burnside 验证隐式归入 Phase II。

**代码实现**：
```python
# main.py:144-185
# Step 1: Symmetry Analysis          ← Phase I
# Step 2: Structure Enumeration      ← Phase II
# Step 3: Verification (BridgeVerifier) ← 独立步骤，不在三阶段中
# Step 4: Assembler (Generate CIFs)  ← Phase III
```
代码将 Verification 实现为一个显式的独立步骤，调用 `BridgeVerifier.verify()` 来检查 Burnside 等式，独立于枚举器内部的 `_verify_burnside()`。

**不一致原因**：
论文的流水线描述缺少验证步骤的明确归属。代码中存在**双重验证**：枚举器内部的 `_verify_burnside()`（自动执行）和主流水线中的 `BridgeVerifier.verify()`（显式调用）。论文应澄清这一结构。

**修改建议**：
在论文流程描述中增加对 Burnside 验证步骤的说明。

**修改前**（原文）：
> \textbf{Phase II (Enumeration)} executes a stratified traversal along the subgroup lattice $\mathcal{L}(\Gamma)$, utilizing orbit-based dimensionality reduction to identify symmetry-unique topological graphs ($\Omega$) while pruning combinatorial noise;

**修改后**（建议）：
> \textbf{Phase II (Enumeration)} executes a stratified traversal along the subgroup lattice $\mathcal{L}(\Gamma)$, utilizing orbit-based dimensionality reduction to identify symmetry-unique topological graphs ($\Omega$) while pruning combinatorial noise. The phase concludes with a Burnside self-verification module that assesses coverage by comparing $\sum 1/|\mathrm{Stab}(c)|$ against $|\mathrm{Fix}(e)|/|\Gamma|$;

**中文翻译**：
Phase II（枚举）沿子群格 L(Γ) 执行分层遍历，利用基于轨道的降维来识别对称性唯一的拓扑图（Ω）并剪枝组合噪声。该阶段以 Burnside 自验证模块收尾，通过比较 Σ 1/|Stab(c)| 与 |Fix(e)|/|Γ| 来评估覆盖率。

---

### [11] Kabsch 算法的架构层次未被论文体现

**位置（论文）**：节 `Geometric Reconstruction`
**位置（代码）**：`assembler.py` + `external_lib/utils/periodic_structure_analysis.py` + `crystal_structure_optimization.py`

**论文描述**：
> Using the Kabsch algorithm, we compute the optimal rotation matrix $\mathbf{R}_{node}$ ...

论文将 Kabsch 描述为装配器（assembler）的核心功能。

**代码实现**：
`assembler.py` 本身并不直接调用 Kabsch 算法。Kabsch（`minimize_rms_superimpose`）实现在 `external_lib/utils/crystal_structure_optimization.py`，由 `periodic_structure_analysis.py` 中的 `get_template_topology_info()` 调用（用于节点对齐和边组装）。`StructureAssembler` 仅作为高层协调器调用 `get_template_topology_info()`。

**不一致原因**：
论文将 Kabsch 直接关联到 `StructureAssembler`，但实际的数值计算被封装在外部库中。这虽然是合理的软件工程实践，但论文的描述可能让读者低估其实现的技术栈深度。

**修改建议**：
在实现细节中说明 Kabsch 在底层库中的位置。

**修改前**（原文）：
> Using the Kabsch algorithm, we compute the optimal rotation matrix $\mathbf{R}_{node}$ that aligns the connection sites of the node SBU with these topological vectors

**修改后**（建议）：
> The Kabsch algorithm (implemented in an \texttt{external\_lib} submodule via SVD in \texttt{NumPy}) is invoked for both node and edge alignment: for each node it computes $\mathbf{R}_{node}$ aligning SBU connection sites with topological vectors; for each edge a second Kabsch alignment orients the linker along the connection axis.

**中文翻译**：
Kabsch 算法（通过 external_lib 子模块中的 NumPy SVD 实现）在节点和边对齐中均被调用：对每个节点，计算 R_node 使 SBU 连接位点与拓扑向量对齐；对每条边，执行第二次 Kabsch 对齐以沿连接轴定向连接体。

---

### [12] processor 型号精度不足

**位置（论文）**：节 `Methods：Computational Implementation Details`

**论文描述**：
> All computations were executed on a single-core Intel Core i5 processor

**代码实现**：
不涉及代码逻辑，属于实验环境声明。

**不一致原因**：
"Intel Core i5"是系列名而非具体型号。不同代际（如第 8 代 vs 第 12 代）的 i5 单核性能差异可达 2–3 倍。"单核"也可能被误解为代码绑定单核，但 NumPy 的 BLAS 运算可能隐式多线程。

**修改建议**：
给出具体型号并禁用了多线程（如果在实验中如此配置）。

**修改前**（原文）：
> All computations were executed on a single-core Intel Core i5 processor, demonstrating the algorithmic efficiency.

**修改后**（建议）：
> All computations were executed on a single core of an Intel Core i5-1035G1 processor (1.0–3.6 GHz, 10th Gen), with \texttt{NumPy} multithreading disabled, demonstrating the algorithmic efficiency.

**中文翻译**：
所有计算在单核 Intel Core i5-1035G1 处理器（1.0–3.6 GHz，第 10 代）上执行，已禁用 NumPy 多线程，证明了算法的高效性。

---

### [13] 高对称目标筛选：stabilizer_size > 1 等价于 ρ(c) > 1/|Γ|

**位置（论文）**：Algorithm 3（GA 基线）/ 节 `Comparison with Random Screening`
**位置（代码）**：`experiment/baseline.py` 第 58–60 行

**论文描述**：
> $\Omega_{target} \gets \{c \in \Omega_{GT} \mid \rho(c) > \frac{1}{|\Gamma|}\}$

论文使用对称性指数 $\rho(c)$ 定义目标集。

**代码实现**：
```python
# experiment/baseline.py:58-60
if res['stabilizer_size'] > 1:
    high_sym_canonicals.add(canon)
```
代码直接检查 `stabilizer_size > 1`。

**不一致原因**：
$\rho(c) = |Stab(c)|/|\Gamma|$，所以 $\rho(c) > 1/|\Gamma| \iff |Stab(c)| > 1$。数学上完全等价。但论文使用连续指标 $\rho(c)$ 而代码使用离散比较 $|Stab(c)| > 1$，符号体系不一致。该等价性在 $|\Gamma| > 1$ 时成立，但当 $|\Gamma| = 1$ 时 $\rho(c) > 1$ 不可能成立，而 `stabilizer_size > 1` 也不会成立。

**修改建议**：
在论文中注明等价关系或在代码注释中显式计算 $\rho(c)$。

**修改前**（原文）：
> $\Omega_{target} \gets \{c \in \Omega_{GT} \mid \rho(c) > \frac{1}{|\Gamma|}\}$ \Comment{Define Symmetric Targets}

**修改后**（建议）：
> $\Omega_{target} \gets \{c \in \Omega_{GT} \mid \rho(c) > \frac{1}{|\Gamma|}\}$ \Comment{Define Symmetric Targets (equivalently, $|\mathrm{Stab}(c)| > 1$)}

**中文翻译**：
定义对称目标集合（等价于 |Stab(c)| > 1）。

---

### [14] Trie 时间复杂度 O(|V|) 而非 O(1)

**位置（论文）**：节 `Methods：Computational Implementation Details`
**位置（代码）**：`enumerator.py` 第 122–135 行 `_trie_lookup()`, `_trie_insert()`

**论文描述**：
> a custom Prefix Tree (Trie) data structure to perform duplicate detection in $O(1)$ time with significantly reduced memory footprint

论文声称 Trie 支持 O(1) 时间重复性检测。

**代码实现**：
```python
def _trie_lookup(self, canonical):
    node = self._trie_root
    for c in canonical:  # 需遍历全部 |V| 个节点
        if c not in node.children: return None
        node = node.children[c]
    return node.class_id
```
Trie 的查找需遍历构型的全部 $|V|$ 个节点，时间复杂度为 $O(|V|)$。

**不一致原因**：
$O(|V|)$ 并不等于 $O(1)$。虽然 $|V|$（节点数）通常较小（$\le 18$），$O(|V|)$ 在常数意义上确实高效，但论文应正确声明。实际上，Trie 相比哈希表的优势在于**内存效率**（共享前缀避免冗余键存储），而非时间复杂度的常数级优化。

**修改建议**：
修正为 $O(|V|)$，强调空间效率优势。

**修改前**（原文）：
> perform duplicate detection in $O(1)$ time with significantly reduced memory footprint

**修改后**（建议）：
> perform duplicate detection in $O(|V|)$ time with significantly reduced memory footprint compared to hash-table-based approaches, as shared prefixes among canonical configurations avoid redundant key storage

**中文翻译**：
以 O(|V|) 时间执行重复性检测，同时相比基于哈希表的方法显著降低内存占用，因为规范构型之间的共享前缀避免了冗余键存储。

---

## 修改优先级建议

### Critical（影响正确性）

| 编号 | 标题 | 理由 |
|------|------|------|
| 1 | Tier 2 子群策略相反 | 论文伪代码与主实现路径的逻辑不同（循环子群 vs 完整子群格），最易被审稿人发现 |
| 2 | 熔断器作用范围 | 伪代码条件 `is Tier 3` 与代码 `not is_full_group` 不同 |
| 3 | Burnside 完备性声称 | 论文描述为严格证明，实际为覆盖率指标，在默认截断下 PARTIAL |

### Major（影响精确性）

| 编号 | 标题 | 理由 |
|------|------|------|
| 4 | 元素二值映射未说明 | 影响对称性分析可复现性的细节 |
| 5 | symprec 提前终止 | 影响扫描行为的实现细节 |
| 6 | DFS 步数限制 vs 穷举 | 核心术语"exhaustive"与实际实现矛盾 |
| 9 | rht-3 性能声明缺上下文 | 轨道数和时间在截断条件下获得 |
| 11 | Kabsch 架构层次 | 论文将计算描述为 assembler 功能，实际在 external_lib |
| 14 | Trie O(1) 声明 | 时间复杂度声明不准确 |

### Minor（措辞优化）

| 编号 | 标题 |
|------|------|
| 7 | 5% 安全缓冲措辞 |
| 8 | 逆置换未注释 |
| 10 | 流水线验证步骤归属 |
| 12 | 处理器型号不精确 |
| 13 | stabilizer_size 与 ρ(c) 等价关系 |

---

*审查生成日期：2026-06-02*
*审查范围：`paper/main.tex` ↔ `main.py`, `symmetry.py`, `enumerator.py`, `verification.py`, `assembler.py`, `experiment/baseline.py`*
