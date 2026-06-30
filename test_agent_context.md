## 项目核心模块

### 核心算法层

- **`SymmetryEngine` (symmetry.py)** — 基于 pymatgen SpacegroupAnalyzer 计算拓扑图的自同构群（automorphism group），包含 lattice 转置、元素清洗（Zn/C 替换）、多精度回退、镜像过滤等容错策略。
- **`MOFEnumerator` (enumerator.py)** — 实现子群格下降（Subgroup Lattice Descent）DFS 算法，对每组轨道（orbit）尝试所有配体颜色分配，通过度数计数约束和正则对称性剪枝进行早期切割；使用 Trie 去重 + 大顶堆维护 Top-K 结果。
- **`BridgeVerifier` (verification.py)** — 使用 Burnside Lemma 进行算法自验证：Sum(1/|Stab(s)|) 应等于 |Fix(e)| / |G|，分"固定组分"和"所有组分"两类计算理论右侧值。
- **`StructureAssembler` (assembler.py)** — 将枚举的着色方案转换为 3D 晶体结构（CIF），包含晶格自适应缩放、拓扑展开（alpha_matrix/Bstar）、X 原子融合、空间群恢复、CIF 后处理（1-indexed 标签、_geom_bond 表）。
- **`StructureValidator` (structure_validator.py)** — 三层硬性校验：碰撞检测（vdW 半径比）、拓扑还原校验（图同构）、对称性守恒校验（空间群偏离 P1 判定）。
- **`enrich_graph_data_with_topology` (main.py)** — 使用真实的晶格平移 sym 向量计算 BFS 路径构建环基（cycle basis），生成 Bstar 和 alpha_matrix 矩阵，通杀扭曲或多重边 MOF。
- **`run_enumeration` (main.py)** — 管线的参数化入口，串联对称性分析 → 枚举 → 验证 → 组装全流程。
- **`Edge` / `Node` (external_lib/)** — 从 CIF 文件解析配体/节点的连接点坐标、unit cell 参数。
- **`get_template_topology_info` / `calculate_X_X_bonds_and_update_charges` / `adjust_placed_edges_to_fit_periodicity` (external_lib/utils/periodic_structure_analysis.py)** — 拓扑图到 3D 坐标的数学展开、X-X 键计算和周期性修正。

### Streamlit 交互层

- **`app.py`** — Streamlit 前端，提供拓扑/配体/边缘选择、3Dmol 可视化、CIF 下载、JSON 导出功能。自动解析节点连通度并动态生成多选菜单。
- **`generate_renders.py`** — 批量生成 3D 渲染图。
- **`plot_energy_landscape.py`** — 能量景观可视化。
- **`visualize_results.py`** — 结果可视化汇总。

### 测试工具层（tests/）

- **`test_mof_pipeline.py`** — ffc 拓扑的管线集成测试（对称性、枚举、验证、组装）。
- **`test_mof_vba.py`** — vba 拓扑的管线集成测试（同上模式）。
- **`run_test.py`** — 交互式参考 CIF 比对测试，含自动修复 Agent 调用。
- **`comparator.py`** — 参考 CIF vs 生成 CIF 的 7 项结构对比（原子数、晶胞参数、体积、空间群、元素组成、原子间距、成键间距）。
- **`name_parser.py`** — CIF 文件名解析（拓扑名、节点列表、edge 名）。
- **`repair_agent.py`** — 基于 Anthropic API 的内嵌自动修复 Agent（5 步循环：确认根因 → 定位代码 → 确认方案 → 执行补丁 → 重测）。

---

## 关键业务规则

### 对称性与自动形态
1. **镜像操作过滤**：`compute_automorphism_group()` 中，旋转矩阵行列式 < 0 的对称操作（含镜像/反演）直接丢弃，仅保留纯旋转置换。
2. **Lattice 转置约定**：pymatgen 使用行向量 × lattice^T，而拓扑 JSON 的 lattice 矩阵需 `.T` 转置后使用，弄反会导致晶格完全错误。
3. **元素清洗**：节点元素含 "Z" 或 "M" 的统一替换为 "Zn"，其余替换为 "C"——这是为了给 SpacegroupAnalyzer 提供合理输入，不影响实际组装。
4. **距离容差 0.8 Å**：对称操作后匹配节点位置的容忍距离，需 `.flatten()` 确保一维数组。
5. **多精度回退**：symprec 依次尝试 [0.01, 0.1, 0.3]，选择返回对称操作最多的分析器。
6. **恒等变换强制包含**：即使 pymatgen 未返回恒等操作，也必须手动添加 `tuple(range(n))`。

### 子群格下降与枚举约束
7. **全群 (G) → 循环子群 (Cyclic) → 恒等群 (Identity) 三级任务**：按轨道数量升序执行，轨道越少（对称性越高）优先。
8. **空间限制 SPACE_LIMIT = 5×10⁶**：非 Identity 任务如果搜索空间 `|atom|^n_orbits` 超过此阈值直接跳过。
9. **步数限制 DFS_STEP_LIMIT = 200,000**：单次 DFS 最多 20 万步（约 1-2 秒），防止 Identity 兜底任务无限运行。
10. **度数统计剪枝**：每个轨道按 `user_selected_ligands[degree]` 选配体；如果某度数节点已超 `topology_requirements[degree]`，立即剪掉该分支。
11. **正则对称性剪枝**：对于每个活跃置换，检查 `g(coloring) < coloring`（字典序），若成立则当前分支被剪掉（已由更小的代表元覆盖）。

### 结构组装约束
12. **晶格自适应缩放**：当配体半径之和 > 节点间实际距离时，缩放晶格 `max_scale * 1.05`；缩放后的键长容差 `bond_tol = 2.5 * scale`。
13. **P1 空间群强制约束**：若 spglib 对称性检测异常（TypeError / Exception），回退到 P1 空间群。
14. **X 原子融合 (Merge)**：节点侧 X 原子（距 edge 侧 X < 0.8 Å）→ 融合到 edge 侧，保留配体侧坐标，移除节点侧副本，且跳过融合后自环和 < 0.8 Å 的短键。
15. **Bridge Bond 1.485 Å 覆写**：Node↔Edge 桥接键的距离强制设为 1.485 Å，键类型设为 'S'。
16. **merge_sites 容差 0.2 Å**：去除融合后的冗余原子，容差必须为 0.2 而非 1.0（否则会误删正常 C-C/C-H 键原子）。
17. **占位键过滤**：距离 = 1.0 Å 的键为 block 连接代码遗留的占位符，CIF 输出时必须跳过。

### 验证规则
18. **Burnside Lemma 验证**：LHS = Σ(1/|Stab|) 必须等于 RHS = |Fix(e)| / |G|，浮点容差 1e-9。
19. **碰撞判定**：d < 0.6 × (vdW_i + vdW_j) 视为碰撞；已键合的原子对不参与碰撞检测。
20. **拓扑同构**：组装后图必须是原始拓扑图的同构图（节点数、边数、度序列、同构性四重校验）。

---

## 外部依赖

| 库 | 用途 | 是否需 Mock |
|----|------|-------------|
| `pymatgen` | Structure/Lattice/ SpacegroupAnalyzer — 结构构建和对称性分析的核心引擎 | **必须 Mock**（对称性分析调用 spglib C 扩展，高耗时且依赖外部库版本） |
| `numpy` | 矩阵运算、线性代数（lattice 变换、距离计算） | 不需要 Mock（纯数值计算） |
| `networkx` | MultiGraph 构建、图同构检验 `nx.is_isomorphic()` | 可选择性 Mock（图同构检验在大图上较慢） |
| `matplotlib` | 图表绘制（实验模块） | 不需要 Mock（仅实验/可视化用） |
| `ase` | 原子模拟环境（实验模块） | 不需要 Mock（非核心管线依赖） |
| `chgnet` | 电荷预测网络（实验模块） | 不需要 Mock（非核心管线依赖） |
| `streamlit` | Web 前端框架 | 不需要 Mock（前端测试用 streamlit 自带 test framework） |
| `py3Dmol` / `stmol` | 3D 分子可视化 | 不需要 Mock（前端仅展示） |
| `anthropic` | Claude API（仅 repair_agent.py） | 不需要 Mock（非核心管线） |

**Mock 优先级说明**：`pymatgen.symmetry.analyzer.SpacegroupAnalyzer` 是最大的性能瓶颈——它调用 spglib 进行对称性检测，单次调用可耗时数秒。所有涉及 `SymmetryEngine` 的单元测试都应 mock 其 `get_symmetry_operations()` 返回值。

---

## 复杂/容易出错的地方

### 对称性分析（symmetry.py）
- **Lattice 转置方向**：`raw_lattice.T` vs `raw_lattice` — 搞反会导致晶格完全错误但不会报错，是最隐蔽的 bug 源。
- **`get_all_distances` 二维返回**：新版本 pymatgen 可能返回 shape (1, N) 的二维数组，`.flatten()` 缺失会导致 `np.where` 始终返回行索引 0（"全是 Node 0" bug）。
- **镜像过滤判断**：`np.linalg.det(rotation_matrix) < 0` 即可靠吗？需要考虑数值误差边界。
- **多精度回退的边界**：如果三个 symprec 都例外失败，返回 `None` → `compute_automorphism_group` 返回仅含恒等变换的列表。

### 子群格下降枚举（enumerator.py）
- **状态回溯完整度**：DFS 中 `current_degree_counts` 和 `current_coloring` 都需正确回溯，遗漏会导致计数紊乱或幽灵着色。当前代码在剪枝 `continue` 分支和正常回溯两处各做了一次回溯，容易不对称。
- **`active_perm_indices` 的传递**：每次递归创建新列表 `next_active_indices`，只有 `g_is_larger == False`（即 g(coloring) >= coloring）的置换才继续活跃。理解错 "canonical representation" 的剪枝逻辑会阻塞合法结果。
- **`top_k_heap` + `max_results` 的交互**：堆中保留的是 stabilizer_size 最大的 K 个，但如果 `max_results` 被用于限制而非 Top-K 语义，行为会异常。
- **`enumerate_subgroups` 的提前退出逻辑**：当前代码 `if max_results and len(self.classes) >= max_results:` 下的退出被注释掉了（只有 `pass`），意味着即使已找够仍会继续跑完所有任务——这是有意为之还是遗留 bug？

### 结构组装（assembler.py）
- **X 原子融合的匹配阈值 0.8 Å**：太小会漏融（产生多余 X），太大会误融（删错原子）。
- **`merge_map` 键重定向**：融合后 `fixed_bonds` 中所有引用被融合原子的键需重定向到保留的 edge 侧原子，遗漏会导致孤立原子或断键。
- **`bridge_bonds` 判定**：通过 `fragment_of` 判断两边分属 node 和 edge 块，如果 `fragment_of` 映射不完整，桥接键判定就会漏检。
- **alpha_matrix / Bstar 的数学正确性**：`enrich_graph_data_with_topology` 中的环基计算依赖于 `omega` 矩阵（边向量）和 `C` 矩阵（环-边关联），sym 解析中的 `int(s_str[0]) - 5` 偏移假设较为脆弱。
- **CIF 标签 0-indexed → 1-indexed 转换**：正则替换 `\b([A-Z][a-z]?)(\d+)\b` → 数字+1，可能错误匹配元素符号。

### 结构校验（structure_validator.py）
- **`bond_exclusions` 构建**：依赖 placed_all 坐标与 structure 坐标精确匹配（< 0.01 Å），若坐标经过 `adjust_placed_edges_to_fit_periodicity` 修改，可能匹配失败导致碰撞检测误报。
- **拓扑还原校验**：节点间通过 edge 块间接连接的推导逻辑（对每个 edge 块的邻居取对称差）假设了 edge 块恰好连接两个 node 块。
- **对称性守恒校验**：在 `_check_symmetry` 中重新调用 `SymmetryEngine` 构建简化图数据，但仅当节点 ID 均为数字时才启用，存在静默跳过的分支。

### 前端交互（app.py）
- **Streamlit 状态管理**：`st.session_state.run_results` 跨 re-run 保持，但不含错误恢复状态；枚举失败后部分日志可能已写入但结果对象不完整。
- **`user_selected_ligands` 构建时机**：依赖 `calculate_topology_requirements` 的返回值动态生成 multiselect，如果拓扑文件缺失或格式异常会导致静默空白。

---

## 已有测试

### 已有测试覆盖

- **`tests/test_mof_pipeline.py`**（ffc 拓扑）
  - `test_symmetry_engine`：验证自同构群非空且含恒等变换
  - `test_mof_enumerator`：验证枚举产出的结构含 coloring / stabilizer_size（限制 5 条结果）
  - `test_verification_logic`：验证 Burnside Lemma LHS=RHS（单配体全穷举）
  - `test_assembler_integration`：集成测试完整 CIF 生成（含 tmp_path）
- **`tests/test_mof_vba.py`**（vba 拓扑）
  - 四个测试用例与 ffc 测试结构一致，但 vba 的 assembler 测试额外运行了枚举来获取着色方案
- **`tests/name_parser.py`** 含内联单测（`_test()`）验证 6 种文件名格式解析
- **`tests/run_test.py`** 交互式对比测试（需人工操作）
- **`tests/repair_agent.py`** 自动修复 Agent（需 ANTHROPIC_API_KEY）

### 明显测试盲区

1. **核心算法层**
   - `enrich_graph_data_with_topology` 的 BFS 环基计算完全没有测试——这是 alpha_matrix 正确性的根本保障
   - `MOFEnumerator._dfs_lattice_pruning` 的剪枝逻辑没有独立单元测试（度数剪枝、正则剪枝的正确性）
   - `MOFEnumerator._get_canonical` / `_compute_stabilizer_and_orbit` 没有独立测试
   - `MOFEnumerator._try_add_result` 的 Top-K 堆管理逻辑（边界：稳等器相等时哪个被替换？）
   - `SymmetryEngine._get_robust_analyzer` 多精度回退逻辑
   - `BridgeVerifier._nCr_multi` 多项式系数计算（边界条件如 count 全 1、count 全 N）

2. **结构组装层**
   - `StructureAssembler._compute_lattice_scale` 晶格缩放计算（边界：无 edge 可用、forced_edge 为空）
   - `StructureAssembler.assemble_and_save` 的 X 原子融合逻辑
   - `_geom_bond` 表的 bridge bond 1.485 Å 覆写逻辑
   - CIF 后处理的 0→1 indexed 标签转换
   - `adjust_placed_edges_to_fit_periodicity` 中 `node_connection_points` 为空时的回退行为

3. **结构校验层**
   - `StructureValidator._check_clashes` 碰撞检测（边界：全 H 结构、超大结构 > 500 原子）
   - `StructureValidator._check_topology` 的 node-edge-node 间接连接推导
   - `StructureValidator._check_symmetry` 的 P1 回退逻辑（边界：拓扑高对称但结构 P1）

4. **边界/异常情况**
   - **空图/单节点拓扑**（仅 1 个节点，恒等群，无需枚举）
   - **无对称性拓扑**（自同构群 = {Identity}，枚举应等价于全排列）
   - **所有原子类型相同**（只有一种配体，枚举结果必唯一）
   - **配体库不匹配拓扑需求**（拓扑需要度-5 节点，但 nodes/ 中无 5c_ 配体）
   - **Graph JSON 文件损坏或缺失关键字段**（lattice/nodes/edges）
   - **CIF 文件不含 X 原子**（节点无连接点）
   - **SpacegroupAnalyzer 全部 symprec 都报异常的退化情况**

5. **跨拓扑泛化**
   - 目前仅测试 ffc 和 vba 两种拓扑，data/ 目录下还有 20+ 种拓扑（cdl, csq, fcu, pcu, sod, tbo 等），从未覆盖
   - 不同拓扑的节点度分布差异大（3-c 到 12-c），未验证枚举器对各度数是否通用

6. **前端测试**
   - Streamlit app.py 完全无测试覆盖（UI 组件渲染、交互流、异常路径）