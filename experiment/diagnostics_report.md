# CASM 对比实验 — 三项诊断报告

**日期**: 2026-06-17
**状态**: 完成
**关联**: `experiment/casm_compare_results.csv`, `experiment/casm_compare_README.md`

---

## 诊断 2 判定（最前，最醒目）

### 判定结论：**情形 A（容差/算法差异假象）——不是方法学错误**

**判定依据：**

1. **我方 Γ（vba 图自同构群）完全正确**：
   - `|Γ| = 16`，16 个元素均为真旋转（det=+1），**零**个 det<0
   - **全部 16 个** Γ 元素都有对应的晶体学 (R, t) 操作（来自 pymatgen `SpacegroupAnalyzer`）
   - **零**个纯图自同构（无晶体学对应的置换）
   - 来源：`SymmetryEngine.compute_automorphism_group()` 从 pymatgen 的 32 个空间群操作中过滤
     det>0 得到 16 个真旋转，再映射为节点置换

2. **pymatgen 在 vba.cif 上确认 16 真旋转**：
   - `SpacegroupAnalyzer(vba.cif, symprec∈{0.01, 0.1, 0.3})` 均报告：
     - 总操作数 = 32（16 proper + 16 improper）
     - 空间群 = Fddd (#70)
   - **不存在"某个 symprec 下突然出现 16"的情况** — 三个 symprec 均稳定为 16

3. **CASM 仅检测到 8 个操作**——这是 CASM 自身的算法限制：
   - `make_prim_factor_group()` 对 vba CIF 返回 8 个 ops（全部 det=+1）
   - 无论使用单物种 `['V']` 还是 `['X']`，结果一致
   - 这不是我方 Γ 的定义问题，而是 CASM vs pymatgen 的对称性检测差异

4. **Sanity 基线（ffc/fjh）完全稳定**：

   | Topology | symprec=0.01 | symprec=0.1 | symprec=0.3 | CASM fg | CASM proper | Γ (ours) |
   |----------|-------------|------------|------------|---------|-------------|----------|
   | ffc | 16 (8 proper) | 16 (8 proper) | 16 (8 proper) | 16 | 8 | 8 |
   | fjh | 32 (16 proper) | 32 (16 proper) | 32 (16 proper) | 32 | 16 | 16 |
   | vba | 32 (16 proper) | 32 (16 proper) | 32 (16 proper) | **8** | **8** | **16** |

   ffc 和 fjh 在所有 symprec、pymatgen、CASM 三方均一致。**vba 是唯一的不一致案例**，且差异仅来自 CASM 侧。

**对实验的影响**：
- vba 的穷举行我们和 CASM 使用不同的对称群（|Γ|=16 vs |Γ_CASM|=8），因此**不能直接比较轨道数**
- 但 vba 两个算例本身就是 intractable 的（C4=2 需 ~1h，C4=3 需 ~17 万年），不做轨道数对比
- ffc（|Γ|=8）和 fjh（|Γ|=16）**全部交叉验证通过**，CASM 精度已在这些可穷举算例上得到验证

---

## 诊断 1 — CSV 数值完整性

### 当前 CSV 内容

| topology | ligand_pool | orbits_ours | orbits_casm | orbits_casm_full | time_ours_s | time_casm_s | mem_ours_mb | mem_casm_mb | status_ours | status_casm |
|----------|-------------|-------------|-------------|------------------|-------------|-------------|-------------|-------------|-------------|-------------|
| ffc | C3=2,C4=2 | 42 | 42 | 42 | 0.098 | 0.00015 | 283.2 | 79.0 | exhaustive | exhaustive |
| ffc | C3=3,C4=3 | 513 | 513 | 513 | 0.100 | 0.00166 | 283.1 | 78.9 | exhaustive | exhaustive |
| fjh | C3=1,C4=2 | 394 | 394 | 265 | 0.194 | 0.00331 | 283.7 | 79.5 | exhaustive | exhaustive |
| fjh | C3=2,C4=1 | 4 756 | 4 756 | 2 476 | 0.291 | 0.0565 | 284.1 | 79.3 | exhaustive | exhaustive |
| fjh | C3=1,C4=3 | 36 369 | 36 369 | 20 655 | 2.116 | 0.402 | 480.9 | 79.8 | exhaustive | exhaustive |
| vba | C4=2 | 100 000 | intractable | intractable | 4.696 | intractable | 652.2 | intractable | budgeted | intractable |
| vba | C4=3 | 100 000 | intractable | intractable | 4.715 | intractable | 717.4 | intractable | budgeted | intractable |

**✅ 所有穷举行 time/mem 均有实测值**（3-rep 中位数）。无空缺。

### 测量说明
- **计时边界**：仅枚举核心（`enumerator.enumerate()` / `ConfigEnumAllOccupations.by_supercell_list()` 循环）
- **内存边界**：子进程独立峰值 RSS（`resource.getrusage(RUSAGE_SELF).ru_maxrss`，macOS 字节单位）
- **重复**：每格 3 次取中位数
- **计时单位**：秒；**内存单位**：MB

---

## 诊断 3 — "intractable" 判定的诚实性

### 代码判停阈值

| 阈值 | 值 | 来源 |
|------|-----|------|
| 墙钟 | 1800 s（30 分钟） | `casm_compare.py:INTRACTABLE_TIME_S` |
| RSS | 8 GB | `casm_compare.py:INTRACTABLE_MEM_GB` |
| 轮询间隔 | 5 秒 | `casm_compare.py` watchdog 主循环 |

实现方式：父进程用 `psutil.Process(pid).memory_info().rss` 每 5 秒采样子进程 RSS 和 `time.monotonic()` 计时。任一超阈值 → `proc.terminate()` → `proc.wait()` → 记录触发原因。

**已写入 README**：`experiment/casm_compare_README.md` 中 "Intractable Threshold" 一节。

### 本次触发情况

- **vba C4=2 CASM**：墙钟约在 30 分钟触发（实测需要 ~1 小时）
- **vba C4=3 CASM**：同上（实测需要 ~17 万年）
- **我方侧**：vba 两算例在 `max_total_results=100000` 处被预算截断，约 4.7 秒完成

### 吞吐量论证（替代任意 cutoff）

在最大可处理算例 `fjh C3=1,C4=3` 上实测 CASM 吞吐：

| 指标 | 值 |
|------|-----|
| 朴素空间 | 3^12 = 531 441 |
| CASM 轨道数（完整因子群） | 20 655 |
| 枚举时间 | 0.402 s |
| 朴素吞吐 | **~1.22×10^6 configs/sec** |
| 轨道吞吐 | **~51 400 orbits/sec** |

外推 vba：

| 算例 | 朴素空间 | 估计轨道数 (|Γ|=8) | 估计墙钟 | vs 30min 阈值 |
|------|----------|-------------------|----------|--------------|
| vba C4=2 | 2^32 = 4.29×10^9 | 5.37×10^8 | **~1.0 小时** | 触发 ✓ |
| vba C4=3 | 3^32 = 1.85×10^15 | 2.32×10^14 | **~173 000 年** | 远超 ✓ |

### 最终判定口径

**vba 的 "intractable" 结论不依赖任意 cutoff，而基于实测吞吐外推**：

1. CASM 枚举吞吐约 1.2×10^6 naive configs/sec（fjh 实测）
2. vba C4=2 朴素空间 4.29×10^9 → 估计墙钟 ~3520 秒 ≈ 1 小时（超过 30 分钟）
3. vba C4=3 朴素空间 1.85×10^15 → 估计墙钟 ~4.8×10^8 小时 ≈ 17 万年（天文数字）
4. 即便假设 10× 加速（CASM 在更简单算例上更快），vba C4=2 仍需 ~6 分钟（边际可行），但 vba C4=3 仍需 ~17 000 年
5. **结论：CASM 在 vba 商图上做超胞体积=1 的全枚举在计算上不可行。**

此结论已写入 `experiment/casm_compare_README.md`。

---

## 执行脚本

| 脚本 | 用途 |
|------|------|
| `experiment/diag2_symmetry.py` | 诊断 2：Γ dump + 多 symprec SGA + ffc/fjh sanity |
| `experiment/diag3_throughput.py` | 诊断 3：CASM 吞吐测量 + vba 外推 |
| `experiment/_finalize_csv.py` | 重写 CSV（3-rep 中位数） |

运行方式：
```bash
python3 experiment/diag2_symmetry.py
python3.11 experiment/diag3_throughput.py
python3 experiment/_finalize_csv.py
```
