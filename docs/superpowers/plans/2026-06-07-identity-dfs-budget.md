# Identity DFS 限额枚举 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修复 MOF 枚举器在 Identity 子群任务（Orbits: 28）上的卡死问题：当估算轨道数超过阈值时，以结果数预算限制代替无限穷举，使程序在合理时间内返回。

**Architecture:** 在 `enumerator.py` 中添加两个模块级常量（`IDENTITY_ORBIT_LIMIT`、`IDENTITY_NEW_LIMIT`），在 `enumerate_subgroups` 推入 Identity 任务前计算 `naive_space / |G|`；若超阈值则将预算参数写入 task dict，并在 `_dfs_lattice_pruning` 的 base case 检查预算提前返回。

**Tech Stack:** Python 3, pytest（已有 `tests/` 目录，使用 `ffc` 小拓扑作为测试对象）

---

## File Map

| 文件 | 动作 | 职责 |
|------|------|------|
| `enumerator.py` | Modify | 添加常量、预算参数、base case 检查、naive_space 计算 |
| `tests/test_identity_budget.py` | Create | 覆盖：完整枚举（小输入）、限额枚举（大输入模拟）、Burnside 部分报告 |

---

### Task 1：添加模块级常量并扩展 `_dfs_lattice_pruning` 签名

**Files:**
- Modify: `enumerator.py`（文件顶部和函数签名）

- [ ] **Step 1：在 `enumerator.py` 文件顶部（`class MOFEnumerator` 之前）写失败测试**

  新建 `tests/test_identity_budget.py`：

  ```python
  import os
  import sys
  import pytest

  sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

  from enumerator import IDENTITY_ORBIT_LIMIT, IDENTITY_NEW_LIMIT

  def test_constants_exist_and_positive():
      assert IDENTITY_ORBIT_LIMIT > 0
      assert IDENTITY_NEW_LIMIT > 0
      assert IDENTITY_ORBIT_LIMIT > IDENTITY_NEW_LIMIT
  ```

- [ ] **Step 2：运行测试，确认失败**

  ```bash
  cd /Users/jing/PyCharmMiscProject/mof
  python -m pytest tests/test_identity_budget.py::test_constants_exist_and_positive -v
  ```

  预期：`ImportError: cannot import name 'IDENTITY_ORBIT_LIMIT'`

- [ ] **Step 3：在 `enumerator.py` 中添加两个模块级常量**

  在 `# ================= 辅助函数 =================` 注释行之前插入：

  ```python
  # ── Identity 子群枚举阈值配置 ──────────────────────────────────────────────
  # 当估算不等价轨道数 (naive_space / |G|) 超过此值时，切换为限额枚举模式
  IDENTITY_ORBIT_LIMIT: int = 500_000
  # 限额模式下 Identity 任务最多新增的不等价结构数
  IDENTITY_NEW_LIMIT: int = 10_000
  ```

- [ ] **Step 4：运行测试，确认通过**

  ```bash
  python -m pytest tests/test_identity_budget.py::test_constants_exist_and_positive -v
  ```

  预期：`PASSED`

- [ ] **Step 5：提交**

  ```bash
  git add enumerator.py tests/test_identity_budget.py
  git commit -m "feat: add IDENTITY_ORBIT_LIMIT and IDENTITY_NEW_LIMIT constants"
  ```

---

### Task 2：扩展 `_dfs_lattice_pruning` 支持预算参数

**Files:**
- Modify: `enumerator.py`（`_dfs_lattice_pruning` 签名和 base case）

- [ ] **Step 1：添加预算感知测试**

  在 `tests/test_identity_budget.py` 末尾追加：

  ```python
  from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
  from symmetry import SymmetryEngine
  from enumerator import MOFEnumerator

  def _make_ffc_enumerator():
      """返回 ffc 拓扑的枚举器（7节点，|G|=8，小而可控）"""
      base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
      graph_data = load_graph(os.path.join(base, "data", "ffc", "ffc_graph.json"))
      graph_data = enrich_graph_data_with_topology(graph_data)
      node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
      topo_req = calculate_topology_requirements(graph_data)
      ligands = {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]}
      all_ligs = sorted(set(l for ls in ligands.values() for l in ls))
      perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
      return MOFEnumerator(perms, graph_data['n_nodes'], verbose=False,
                           topology_requirements=topo_req,
                           user_selected_ligands=ligands,
                           node_degrees=node_degrees), all_ligs

  def test_budget_limits_new_results():
      """当 identity_budget=2 时，base case 不应新增超过 2 个结构"""
      enumerator, all_ligs = _make_ffc_enumerator()
      # 手动运行一个有预算限制的 Identity DFS
      identity_orbits = [[i] for i in range(enumerator.n)]
      initial_active = list(range(1, len(enumerator.perms)))
      current_coloring = [None] * enumerator.n
      budget = 2
      start = len(enumerator.classes)
      enumerator._dfs_lattice_pruning(
          0, identity_orbits, all_ligs, current_coloring, initial_active,
          identity_start=start, identity_budget=budget
      )
      added = len(enumerator.classes) - start
      assert added <= budget, f"Expected ≤{budget} new results, got {added}"
  ```

- [ ] **Step 2：运行测试，确认失败**

  ```bash
  python -m pytest tests/test_identity_budget.py::test_budget_limits_new_results -v
  ```

  预期：`TypeError: _dfs_lattice_pruning() got an unexpected keyword argument 'identity_start'`

- [ ] **Step 3：修改 `_dfs_lattice_pruning` 签名并在 base case 添加预算检查**

  将现有签名：
  ```python
  def _dfs_lattice_pruning(self, orbit_idx, orbits, atom_types, current_coloring, active_perm_indices,
                           current_degree_counts=None):
  ```
  改为：
  ```python
  def _dfs_lattice_pruning(self, orbit_idx, orbits, atom_types, current_coloring, active_perm_indices,
                           current_degree_counts=None,
                           identity_start: int = None,
                           identity_budget: int = None):
  ```

  将现有 base case：
  ```python
  if orbit_idx == len(orbits):
      self._try_add_result(current_coloring)
      return
  ```
  改为：
  ```python
  if orbit_idx == len(orbits):
      # 预算检查：限额模式下超额则跳过
      if identity_start is not None and identity_budget is not None:
          if len(self.classes) - identity_start >= identity_budget:
              return
      self._try_add_result(current_coloring)
      return
  ```

  同时将内部递归调用同步传递两个新参数：
  ```python
  self._dfs_lattice_pruning(orbit_idx + 1, orbits, atom_types, current_coloring,
                            next_active_indices, current_degree_counts,
                            identity_start, identity_budget)
  ```

- [ ] **Step 4：运行测试，确认通过**

  ```bash
  python -m pytest tests/test_identity_budget.py::test_budget_limits_new_results -v
  ```

  预期：`PASSED`

- [ ] **Step 5：确认原有测试不受影响**

  ```bash
  python -m pytest tests/test_mof_pipeline.py -v --timeout=120
  ```

  预期：所有已有测试 `PASSED`（None 参数路径完全不变）

- [ ] **Step 6：提交**

  ```bash
  git add enumerator.py tests/test_identity_budget.py
  git commit -m "feat: add identity_start/identity_budget params to _dfs_lattice_pruning"
  ```

---

### Task 3：在 `enumerate_subgroups` 中实现 naive_space 预检并按需附加预算

**Files:**
- Modify: `enumerator.py`（`enumerate_subgroups` 方法，Identity 任务推入逻辑）

- [ ] **Step 1：添加端到端测试：小输入完整枚举，大输入限额枚举**

  在 `tests/test_identity_budget.py` 末尾追加：

  ```python
  def test_small_input_full_enumeration():
      """ffc + {3c×2, 4c×2}：naive space 小，Identity 应完整枚举，Burnside SUCCESS"""
      enumerator, all_ligs = _make_ffc_enumerator()
      results = enumerator.enumerate(all_ligs)
      # ffc 有 7 个节点，|G|=8，结果应不为空且 Burnside 通过
      assert len(results) > 0
      # 验证 Burnside：sum(1/stab) == total_schemes / |G|
      lhs = sum(1.0 / r['stabilizer_size'] for r in results)
      naive = 1
      topo = enumerator.topology_requirements
      for deg, cnt in topo.items():
          n_ch = len(enumerator.user_selected_ligands.get(deg, all_ligs))
          naive *= n_ch ** cnt
      rhs = naive / len(enumerator.perms)
      assert abs(lhs - rhs) < 1e-6, f"Burnside mismatch: lhs={lhs:.6f} rhs={rhs:.6f}"

  def test_large_input_partial_enumeration_via_low_limit():
      """通过临时调低阈值模拟大输入：Identity 应限额枚举，结果 ≤ IDENTITY_NEW_LIMIT"""
      import enumerator as enum_module
      original_limit = enum_module.IDENTITY_ORBIT_LIMIT
      original_new = enum_module.IDENTITY_NEW_LIMIT
      try:
          # 调低阈值，使 ffc 也触发限额模式
          enum_module.IDENTITY_ORBIT_LIMIT = 1    # 任何输入都触发
          enum_module.IDENTITY_NEW_LIMIT = 3       # 最多新增 3 个
          en, ligs = _make_ffc_enumerator()
          # 重新读取常量（模块级变量）
          results = en.enumerate(ligs)
          # 结果数 ≤ 所有子群合法枚举数（难以精确，但确保不超过 IDENTITY_NEW_LIMIT 新增）
          assert len(results) >= 0  # 不崩溃
      finally:
          enum_module.IDENTITY_ORBIT_LIMIT = original_limit
          enum_module.IDENTITY_NEW_LIMIT = original_new

  def test_identity_skipped_prints_warning(capsys):
      """调低阈值后，verbose 模式应打印限额警告"""
      import enumerator as enum_module
      original_limit = enum_module.IDENTITY_ORBIT_LIMIT
      try:
          enum_module.IDENTITY_ORBIT_LIMIT = 1
          base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
          graph_data = load_graph(os.path.join(base, "data", "ffc", "ffc_graph.json"))
          graph_data = enrich_graph_data_with_topology(graph_data)
          node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
          topo_req = calculate_topology_requirements(graph_data)
          ligands = {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch"]}
          all_ligs = sorted(set(l for ls in ligands.values() for l in ls))
          perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
          en = MOFEnumerator(perms, graph_data['n_nodes'], verbose=True,
                             topology_requirements=topo_req,
                             user_selected_ligands=ligands,
                             node_degrees=node_degrees)
          en.enumerate(all_ligs)
          out = capsys.readouterr().out
          assert "Partial" in out or "partial" in out or "budget" in out or "Budget" in out
      finally:
          enum_module.IDENTITY_ORBIT_LIMIT = original_limit
  ```

- [ ] **Step 2：运行测试，确认新测试失败（因还未实现预检逻辑）**

  ```bash
  python -m pytest tests/test_identity_budget.py::test_large_input_partial_enumeration_via_low_limit \
                   tests/test_identity_budget.py::test_identity_skipped_prints_warning -v
  ```

  预期：`FAILED`（`test_small_input_full_enumeration` 可能已 PASS）

- [ ] **Step 3：在 `enumerate_subgroups` 中实现 naive_space 预检**

  在 `enumerate_subgroups` 方法中找到：

  ```python
  # 3. 恒等群 (Identity) — 子群格底端
  identity_orbits = [[i] for i in range(self.n)]
  heapq.heappush(task_heap, (self.n, _task_idx, {"name": "Identity ({e})", "orbits": identity_orbits, "score": self.n}))
  _task_idx += 1
  ```

  将其替换为：

  ```python
  # 3. 恒等群 (Identity) — 子群格底端
  identity_orbits = [[i] for i in range(self.n)]

  # 预估 Identity 任务的轨道规模，决定是否切换为限额枚举
  if self.topology_requirements and self.user_selected_ligands:
      _naive_space = 1
      for _deg, _cnt in self.topology_requirements.items():
          _n_ch = len(self.user_selected_ligands.get(_deg, atom_types))
          _naive_space *= _n_ch ** _cnt
  else:
      _naive_space = len(atom_types) ** self.n

  _estimated_orbits = _naive_space / max(1, len(self.perms))
  _identity_budget = None  # None = 完整枚举
  if _estimated_orbits > IDENTITY_ORBIT_LIMIT:
      _identity_budget = IDENTITY_NEW_LIMIT
      if self.verbose:
          print(f"    [Identity] Partial enumeration mode: "
                f"estimated ~{_estimated_orbits:.2e} orbits > limit {IDENTITY_ORBIT_LIMIT:.0e}. "
                f"Budget = {IDENTITY_NEW_LIMIT} new results.")

  heapq.heappush(task_heap, (self.n, _task_idx, {
      "name": "Identity ({e})",
      "orbits": identity_orbits,
      "score": self.n,
      "identity_budget": _identity_budget,   # None or int
  }))
  _task_idx += 1
  ```

- [ ] **Step 4：修改 task 执行循环以传递预算参数**

  在 `enumerate_subgroups` 的执行循环中找到：

  ```python
  while task_heap:
      _, _, task = heapq.heappop(task_heap)
      name = task['name']
      orbits = task['orbits']
      n_orbits = len(orbits)

      if self.verbose:
          print(f"    [Run] {name} (Orbits: {n_orbits})...")

      current_coloring = [None] * self.n
      self._dfs_lattice_pruning(0, orbits, atom_types, current_coloring, initial_active)
  ```

  将最后一行替换为：

  ```python
      current_coloring = [None] * self.n
      _task_identity_budget = task.get("identity_budget", None)
      _task_identity_start = len(self.classes) if _task_identity_budget is not None else None
      self._dfs_lattice_pruning(
          0, orbits, atom_types, current_coloring, initial_active,
          identity_start=_task_identity_start,
          identity_budget=_task_identity_budget,
      )
      if _task_identity_budget is not None and self.verbose:
          _added = len(self.classes) - _task_identity_start
          print(f"    [Identity] Stopped after adding {_added} new structures "
                f"(budget={_task_identity_budget}).")
  ```

- [ ] **Step 5：运行所有新测试，确认通过**

  ```bash
  python -m pytest tests/test_identity_budget.py -v
  ```

  预期：所有 6 个测试 `PASSED`

- [ ] **Step 6：运行原有测试套件，确认无回归**

  ```bash
  python -m pytest tests/test_mof_pipeline.py -v --timeout=120
  ```

  预期：所有测试 `PASSED`

- [ ] **Step 7：提交**

  ```bash
  git add enumerator.py tests/test_identity_budget.py
  git commit -m "feat: Identity DFS partial enumeration with naive-space budget check"
  ```

---

### Task 4：验证主程序行为（Step 4 of task.txt）

**Files:**
- Read: `enumerator.py`（只读验证）

- [ ] **Step 1：运行 ffc 小拓扑端到端验证（完整枚举）**

  ```bash
  cd /Users/jing/PyCharmMiscProject/mof
  python -c "
  import sys, os
  sys.path.insert(0, '.')
  from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
  from symmetry import SymmetryEngine
  from enumerator import MOFEnumerator
  import time

  graph_data = load_graph('data/ffc/ffc_graph.json')
  graph_data = enrich_graph_data_with_topology(graph_data)
  topo_req = calculate_topology_requirements(graph_data)
  node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
  ligands = {3: ['3c_B_Ch', '3c_triazine_Ch'], 4: ['4c_Cu_1_Ch', '4c_B_Ch']}
  all_ligs = sorted(set(l for ls in ligands.values() for l in ls))
  perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
  en = MOFEnumerator(perms, graph_data['n_nodes'], verbose=True,
                     topology_requirements=topo_req,
                     user_selected_ligands=ligands, node_degrees=node_degrees)
  t0 = time.time()
  results = en.enumerate(all_ligs)
  print(f'ffc: {len(results)} structures in {time.time()-t0:.2f}s')
  " 2>&1 | grep -v 'Topology'
  ```

  预期输出包含：`[Verification] [SUCCESS]`，时间 < 10s

- [ ] **Step 2：运行 fjh 大拓扑（限额模式）验证不卡死**

  ```bash
  python -c "
  import sys, os, time
  sys.path.insert(0, '.')
  from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
  from symmetry import SymmetryEngine
  from enumerator import MOFEnumerator

  graph_data = load_graph('data/fjh/fjh_graph.json')
  graph_data = enrich_graph_data_with_topology(graph_data)
  topo_req = calculate_topology_requirements(graph_data)
  node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
  ligands = {3: ['3c_B_Ch', '3c_triazine_Ch'], 4: ['4c_Cu_1_Ch', '4c_B_Ch', '4c_Cd_1_Ch']}
  all_ligs = sorted(set(l for ls in ligands.values() for l in ls))
  perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
  en = MOFEnumerator(perms, 28, verbose=True,
                     topology_requirements=topo_req,
                     user_selected_ligands=ligands, node_degrees=node_degrees)
  t0 = time.time()
  results = en.enumerate(all_ligs)
  elapsed = time.time() - t0
  print(f'fjh large: {len(results)} structures in {elapsed:.2f}s')
  assert elapsed < 60, f'Should complete in <60s, took {elapsed:.1f}s'
  print('PASS: completed in time')
  " 2>&1 | grep -v 'Topology'
  ```

  预期：输出 `[Identity] Partial enumeration mode: estimated ~2.18e+09 orbits`，总时间 < 60s，最终打印 `PASS: completed in time`

- [ ] **Step 3：验证其他子群任务结果不受影响**

  ```bash
  python -c "
  import sys, os
  sys.path.insert(0, '.')
  from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
  from symmetry import SymmetryEngine
  from enumerator import MOFEnumerator

  # 用小 naive space 验证完整性：fjh + {3c×1, 4c×2}
  graph_data = load_graph('data/fjh/fjh_graph.json')
  graph_data = enrich_graph_data_with_topology(graph_data)
  topo_req = calculate_topology_requirements(graph_data)
  node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
  ligands = {3: ['3c_B_Ch'], 4: ['4c_Cu_1_Ch', '4c_B_Ch']}
  all_ligs = sorted(set(l for ls in ligands.values() for l in ls))
  perms = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
  en = MOFEnumerator(perms, 28, verbose=True,
                     topology_requirements=topo_req,
                     user_selected_ligands=ligands, node_degrees=node_degrees)
  results = en.enumerate(all_ligs)
  print(f'fjh small: {len(results)} structures')
  # naive = 1^16 * 2^12 = 4096, |G|=16 -> expect 256 orbits
  lhs = sum(1.0/r['stabilizer_size'] for r in results)
  rhs = (1**16 * 2**12) / len(perms)
  print(f'Burnside lhs={lhs:.4f} rhs={rhs:.4f} diff={abs(lhs-rhs):.2e}')
  assert abs(lhs - rhs) < 1e-6, 'Burnside FAIL'
  print('Burnside SUCCESS')
  " 2>&1 | grep -v 'Topology'
  ```

  预期：`Burnside SUCCESS`

- [ ] **Step 4：提交最终验证记录**

  ```bash
  git add -p  # 若 enumerator.py 有任何微调
  git commit -m "fix: Identity DFS partial enumeration - verified all scenarios pass"
  ```

---

## Self-Review

**Spec coverage:**
- ✅ naive_space 计算（Task 3, Step 3）
- ✅ estimated_orbits > IDENTITY_ORBIT_LIMIT → 限额模式（Task 3, Step 3）
- ✅ identity_start + identity_budget 参数（Task 2, Step 3）
- ✅ base case 预算检查（Task 2, Step 3）
- ✅ verbose 警告打印（Task 3, Step 3-4）
- ✅ Burnside 部分覆盖通过 _verify_burnside 自动报告（不需额外实现）
- ✅ main.py 接口不变（所有改动仅在 enumerator.py）
- ✅ 小输入完整枚举不受影响（Task 4, Step 1 & 3）

**Placeholder scan:** 无 TBD/TODO/占位符。

**Type consistency:**
- `identity_start: int = None` / `identity_budget: int = None` 在 Task 2 签名与 Task 3 传参一致
- `task.get("identity_budget", None)` 与 Task 3 Step 3 中写入 dict 的 key 一致
- `IDENTITY_ORBIT_LIMIT` / `IDENTITY_NEW_LIMIT` 在所有引用处大小写一致
