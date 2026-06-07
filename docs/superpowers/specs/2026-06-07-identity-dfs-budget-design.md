---
title: Identity 子群 DFS 限额部分枚举修复
date: 2026-06-07
status: approved
---

## 问题

MOF 枚举程序在执行子群格下降（`enumerate_subgroups`）时，最后的 Identity 子群任务
（每个节点各自为一个 orbit，共 28 个 orbit）对于大输入卡死无返回。

**根因（情形A：搜索树爆炸）：**

- fjh 拓扑：28 节点（16 个 3c + 12 个 4c），|G| = 16
- 当前配体：{3c: 2 种, 4c: 3 种}
- 朴素搜索空间：2¹⁶ × 3¹² ≈ 348 亿
- Burnside 理论轨道数：≈ 21.8 亿
- Python 遍历速度：< 100 万次/秒 → 需数小时，程序看似卡死

较小输入（如 3c×1, 4c×3 → naive ≈ 531K → 约 3.3 万轨道）可正常完成。

## 设计

### 策略：估算轨道数，超阈值时切换为限额枚举

在 `enumerate_subgroups` 中推入 Identity 任务前：

1. 计算 `naive_space`（复用 `_verify_burnside` 中相同公式）
2. 计算 `estimated_orbits = naive_space / |G|`
3. 比较阈值 `IDENTITY_ORBIT_LIMIT = 500_000`

**分支：**
- `estimated_orbits ≤ IDENTITY_ORBIT_LIMIT`：正常推入 Identity 任务，完整枚举（不变）
- `estimated_orbits > IDENTITY_ORBIT_LIMIT`：推入 Identity 任务，但附带结果数预算：
  - 记录 `identity_start_count = len(self.classes)`
  - `identity_budget = IDENTITY_NEW_LIMIT = 10_000`
  - `_dfs_lattice_pruning` 的 base case 中检查预算，达到则立即返回

### 常量（模块级）

| 常量 | 默认值 | 说明 |
|------|--------|------|
| `IDENTITY_ORBIT_LIMIT` | `500_000` | 超过时切换为限额模式 |
| `IDENTITY_NEW_LIMIT` | `10_000` | 限额模式下最多新增结构数 |

### 接口兼容性

- **不修改** `main.py` 调用接口
- `_dfs_lattice_pruning` 新增两个关键字参数（均默认 `None`）：
  - `identity_start: Optional[int] = None`
  - `identity_budget: Optional[int] = None`
- 当两者均为 `None` 时，行为与原代码完全相同

### 输出变化

限额模式下额外打印：
```
    [Identity] Partial enumeration: estimated ~2.18e+09 orbits, budget=10000 new results.
    [Identity] Stopped after adding N new structures.
```

`_verify_burnside` 自动报告 PARTIAL 覆盖率，不需额外修改。

## 验证标准

1. fjh + {3c×2, 4c×3}：Identity 任务在 <60s 内完成（限额枚举后停止）
2. fjh + {3c×1, 4c×3}：Identity 任务完整枚举，Burnside 验证 SUCCESS
3. ffc（小拓扑）：所有配置均完整枚举，结果数不变
4. main.py 调用接口无变化
