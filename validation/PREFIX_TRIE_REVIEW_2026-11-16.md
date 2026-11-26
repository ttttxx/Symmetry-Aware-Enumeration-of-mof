# 已完成前缀 Trie：实现、核验与效果

## 本轮修改

保留子群 DFS 和当前完整代表去重。新增独立的搜索覆盖 Trie，Python 与
C++ DFS 均可用；`CasmInspiredMOFEnumerator` 原有的完整代表 Trie 也继续保留，
`HashStorageMOFEnumerator`/混合后端的完整字节键哈希未改动。

一次 DFS 中每个前缀通常只进入一次，单纯记录“见过”无法省掉这次搜索。
新索引记录的是：某个任务约束下，这个前缀对应的规范搜索已经全部结束。
完整任务已有整域排除；本轮增加的能力主要是复用**中断任务中已完成的部分**。

实现位于 `search_prefix_trie.py`、`enumerator_casm_inspired.py`、
`enumerator_native.py` 和 `native/enumeration_core.cpp`：

- DFS 保存各层已完成选项数；中断时把当前路径的已完成兄弟分支压缩进 Trie，
  未完成路径保持开放。不在每个搜索节点上分配 Trie 节点。
- 后续任务查询时先确认允许标签域包含关系，并证明当前分支的所有补全均满足
  来源任务的轨道等色约束，然后匹配已关闭前缀。
- 只在相关前缀和等式首次可判定的深度查询。失败的条件在其后代中不会恢复，
  因此同一证书无需逐层重复查询。
- 默认最多缓存 50,000 个节点；容量不足时少存覆盖证书，继续正常搜索。
  该限制针对 Trie 节点数，不是进程内存上限。每次枚举重新建缓存；原生上下文
  销毁时释放 C++ 节点。
- 完整任务继续走原有整域排除，不建前缀节点。全局结果上限停止后通常没有后续
  任务可复用，不能指望该情形自动提速。精确 rho 层筛选尚未集成。

CLI/API 已有显式开关，默认关闭；原默认 legacy 后端保持不变：

```sh
.venv/bin/python tools/build_enumeration_core.py
.venv/bin/python main.py --enumeration-backend cpp --prefix-trie
# Python 对照：将 cpp 换为 python
```

`run_enumeration(..., enumeration_backend="cpp", prefix_trie=True)` 或构造器
`HybridMOFEnumerator(..., backend="cpp", prefix_trie=True)` 同样可用。
节点上限为 `--prefix-trie-max-nodes` / `prefix_trie_max_nodes`。
共享库 ABI 从 1 升为 2，旧库须重建；auto 缺库/版本不符时走 Python 回退。

## 正确性与构建

本机 Clang C++17 构建及 `-Wall -Wextra -Werror -fsyntax-only` 通过。
相关测试合计 **144 passed**，19 条既有 spglib 弃用提示。

新增测试遍历二元四节点任务的全部 0–8 截断预算，分别用 Python/C++ 和
0、2、50,000 节点容量验证；涵盖开放前缀、不同任务划分、标签域不兼容、
单标签常量、变量顺序变化、多次非嵌套任务截断、过期任务、重复枚举、输出上限
和恒等任务预算。独立穷举轨道集合及整数轨道覆盖通过；完整结果列表和稳定子/
轨道大小保持一致。另验证保留完整代表 Trie 时也能使用新前缀缓存。

```sh
.venv/bin/python -m pytest -q tests/test_search_prefix_trie.py \
  tests/test_native_enumerator.py tests/test_enumerator_storage.py \
  tests/test_counter_enumerator.py tests/test_subgroup_task_dedup.py \
  tests/test_identity_budget.py tests/test_casm_scale_benchmark.py
```

固定秒数超时会因运行速度变化返回不同数量结果；本轮未声称固定时间结果集合
相同。截断安全性用确定的新结果预算逐位置核验，过期任务另行测试。

## 同存储 on/off 实验

入口 `validation/benchmark_prefix_trie_20261116.py`；原始结果、计数器、源码及
原生库哈希见 `PREFIX_TRIE_BENCHMARK_2026-11-16.json`。本轮为项目 Python
3.13 环境，同进程预热一轮后轮换五次，记录 `enumerate` 端到端中位数；图/群准备、
GC 和输出核验在计时外，不分析 RSS。7 组 × 4 版本 × 5 次 = 140 次正式计时，
另 28 次预热；每次完整输出顺序及元数据相同，整数轨道覆盖通过。
此口径与之前 CASM 独立进程/Python 3.11 实验不同，不能跨实验直接计算加速比。

人工任务使用仅交换前两节点的群，先搜索其固定点任务，在首个标签对应的整个
子树结束后截断，再搜索恒等任务。原始组合为 2^n，完整代表为 3×2^(n−2)。
该设计让预期能复用的区域和剪枝位置均可确定。

| 情形 | Python 关闭/开启 (s) | C++ 关闭/开启 (s) | DFS 进入节点数 关闭→开启 | 重复叶子 关闭→开启 |
|---|---:|---:|---:|---:|
| 人工 n=12 | 0.007901 / 0.006544 | 0.002873 / 0.002746 | 8,192→6,146 | 1,024→0 |
| 人工 n=16 | 0.139116 / 0.117319 | 0.055966 / 0.052860 | 131,072→98,306 | 16,384→0 |
| 人工 n=18 | 0.804428 / 0.703016 | 0.478972 / 0.461282 | 524,288→393,218 | 65,536→0 |
| ffc 3/3，正常穷尽 | 0.003002 / 0.003093 | 0.002258 / 0.002348 | 1,048→1,048 | 0→0 |
| ffc 3/3，首任务预算 3 | 0.002985 / 0.003531 | 0.002302 / 0.002304 | 1,054→1,051 | 3→0 |
| fjh 1/3，正常穷尽 | 0.560757 / 0.590170 | 0.098012 / 0.096842 | 620,936→620,936 | 0→0 |
| fjh 1/3，首任务预算 3 | 0.560392 / 0.630484 | 0.097652 / 0.097893 | 620,948→620,939 | 3→0 |

“DFS 进入节点数”不包括递归前即被群作用比较拒绝的选项。人工三例均只建两个
Trie 节点，查三次、命中一次，命中位置在完整叶子之前，省掉约 25% 的 DFS
进入次数。端到端 Python 约快 1.14–1.21 倍，C++ 约快 1.04–1.06 倍；C++ 中
输出物化等共同成本占比更高，节点降幅不能直接当作耗时降幅。

正常穷尽的两个真实算例没有前缀节点或新命中：已完整完成的任务本来就受整域
排除保护。首任务预算 3 时两例均一次提前命中，但只省去极少搜索，Python
查询成本反而更高，C++ 总耗时基本不变。因此当前证据支持“可复用较大的已完成
子树时有收益”，不支持“启用 Trie 就会普遍加速”。暂保留可选开关。

后续若继续优化，应先统计真实限时/限量调度里能复用区域的大小及查询成本，
考虑将来源任务的等式失败信息和 Trie 游标沿 DFS 传播，避免约束较晚确定时
重复查询。完整代表存储的进一步改动按作者要求暂缓。

## 历史实验溯源

修改前已将 CASM 梯度实验 manifest 中 13 个测量源码/原生库文件按原哈希归档
到 `validation/snapshots/casm_scale_20261116`。原始 CASM 数值和 manifest 未改，
其时间属于归档的 ABI-1 版本，不代表本轮前缀 Trie 的时间。
分析脚本新增显式来源目录和只核验模式，历史 180 次结果核验通过：

```sh
.venv/bin/python validation/analyze_casm_scale_20261116.py \
  --source-root validation/snapshots/casm_scale_20261116 --check-only
```

本轮未修改论文、正式性能 CSV、CIF/DFT 流程或精确 rho 配额实现。
