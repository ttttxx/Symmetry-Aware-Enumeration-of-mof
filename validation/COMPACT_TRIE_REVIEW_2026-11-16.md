# C++ 紧凑数组 Trie：完整代表去重

## 实现与统一算法

按作者要求将 C++ 完整代表索引切为紧凑数组 Trie，保留同核哈希模式作为对照。
统一的方法定义、状态、完整伪代码、正确性和复杂度边界见
`docs/subgroup_dfs_trie_algorithm.md`，按同一条流程描述：

**子群固定点任务 → 规范前缀 DFS → 已完成区域排除 → 完整代表 Trie 插入 →
稳定子/轨道元数据 → 排序输出。**

完整代表 Trie 和已完成前缀缓存在该流程中各自记录不同事实。完整键出现过
支持精确查重，关闭前缀支持排除整个搜索区域，后者还要求来源任务约束相容。
不将这两种终点混为同一个标记。

`native/enumeration_core.cpp` 新增 `CompactTrie`：

- 路径压缩，连续 `vector<Node>`，每节点 16 字节。
- 四个 32 位字段：样本记录下标、深度、首孩子、下一兄弟。没有每节点 map，
  没有每节点独立堆分配，兄弟按无符号首字节排序。
- 边标签引用原生结果键的区间，不另存标签字符串；插入时在首次差异处分叉。
  所有键固定长度，完整键叶子判重，0 和 128–255 标签均精确支持。
- Trie 模式不往 C++ `unordered_set` 插键，也不构建 Python 完整键哈希副本。
  结束前导出连续节点数组，释放原生上下文后仍可基于结果记录查询。
- 数组和记录下标扩容后保持有效；32 位范围在增加前检查。节点数至多 $2r+1$，
  实际 vector 容量和公共结果键空间另计。

Python `classes` 仍按稳定的插入下标保存记录，返回列表的 `id` 在稳定子排序
后重编号。Trie 返回与既有索引相同的插入下标，不能把它误当作排序后公开 ID。

混合后端默认 `native_store="trie"`，哈希可用 `native_store="hash"` 恢复。
`main.py` 增加 `--native-store`，benchmark JSON 记录请求、实际存储及索引统计。
原生 ABI 升为 3；主 CLI 未指定后端时仍是 legacy，Python 回退保持原有字节/
元组哈希。前缀覆盖仍需 `--prefix-trie` 显式开启。

```sh
.venv/bin/python tools/build_enumeration_core.py
.venv/bin/python main.py --enumeration-backend cpp --native-store trie --prefix-trie
```

## 核验

- 项目环境相关测试 **158 passed**，19 条既有 spglib 弃用提示。覆盖完整结果、
  稳定子、轨道大小、排序后 ID、数量预算、截断前缀复用、缺库回退和索引重置。
- C++ 结构测试使用独立哈希 oracle，键长 1/2/16/512，各六种插入次序，
  每序列 2,000 次主插入并再次插入查重，另随机查询不存在/存在键。检查长公共
  路径、拆边、兄弟顺序、全部节点可达、无环、节点数界及字节 0–255。
  Clang 严格警告、AddressSanitizer 和 UndefinedBehaviorSanitizer 通过。
- Python 查询测试验证解码后使用导出的节点数组、最大 256 标签的全部二元
  指派、缺失键、错误长度、原实例更换标签后的缓存重置，Trie 模式完整键哈希
  条目数为零。正常 DFS 和比较计数器均与原索引输出一致。

```sh
.venv/bin/python -m pytest -q tests/test_native_compact_index.py \
  tests/test_native_enumerator.py tests/test_search_prefix_trie.py \
  tests/test_enumerator_storage.py tests/test_counter_enumerator.py \
  tests/test_subgroup_task_dedup.py tests/test_identity_budget.py \
  tests/test_casm_scale_benchmark.py

clang++ -std=c++17 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Wall -Wextra -Werror \
  tests/native_compact_trie.cpp -o /tmp/mof_compact_trie_test
/tmp/mof_compact_trie_test
```

## 同 C++ DFS、同群、同输出对照

入口 `validation/benchmark_compact_trie_20261116.py`；原始 36 次记录、摘要及
源码/原生库哈希见 `COMPACT_TRIE_BENCHMARK_2026-11-16.json`。六例各三次独立
进程，每轮交替 hash/trie 顺序。项目 Python 3.13.12、C++17 -O3，单线程环境，
强制原生 DFS。两边都开启前缀覆盖；普通穷尽例没有前缀命中，人工例能复用中断
任务的已完成大子树。

图/群分析、SymPy 导入及独立 Burnside 计数在计时外，`enumerate` 的任务准备、
搜索、导出、结果物化、保留索引和排序在计时内。峰值 RSS 在输出摘要计算前
读取，为完整进程高水位，包含公共依赖和结果，不作为纯索引内存。

每例六次运行的完整排序结果、公开 ID 和元数据摘要一致，搜索节点与叶访问数
一致，独立 Burnside 计数和整数轨道覆盖通过。人工任务的首个任务中断、后续
任务完成全空间，最终仍有完整覆盖。

以下均为三次中位数，时间单位秒，RSS 单位 MiB：

| 算例 | 代表数 | 端到端 hash→Trie | C++ 核 hash→Trie | 峰值 RSS hash→Trie |
|---|---:|---:|---:|---:|
| ffc 3/3 | 513 | 0.004218→0.004293 | 0.000133→0.000130 | 131.8→131.9 |
| ffc 5/8 | 58,080 | 0.067244→0.063518 | 0.007367→0.006493 | 169.9→160.8 |
| fjh 1/3 | 36,369 | 0.114647→0.112305 | 0.024833→0.024162 | 164.7→158.1 |
| fjh 1/4 | 1,078,816 | 2.686834→2.118668 | 0.826759→0.636847 | 1123.2→937.6 |
| lil 4/3 | 1,328,832 | 1.995582→1.767814 | 0.378739→0.158133 | 1011.6→871.7 |
| 人工 swap18，部分任务中断 | 196,608 | 0.290119→0.247356 | 0.037753→0.039924 | 306.4→281.6 |

本轮两个百万结果例支持采用紧凑 Trie：fjh 端到端约快 1.27 倍、峰值 RSS 约减
16.5%，lil 约快 1.13 倍、RSS 约减 13.8%。两者 C++ 核本身也更快；端到端收益
还包含取消 Python 哈希副本的效果。人工例原生核略慢、端到端更快，说明二者
应分开解释。小例差异接近固定开销与计时波动，不主张所有输入都提速。

lil Trie 三次端到端为 1.547/1.899/1.768 s，存在明显运行波动，当前报告采用
中位数，不据此建立普遍加速界。此前 Python 对象 Trie 的回退不能直接推到
C++ 压缩数组 Trie；本轮也不以旧 Python Trie 数字计算新实现的加速比。

最大 fjh 有 1,541,373 个索引节点，逻辑节点约 23.52 MiB，vector 保留容量
32 MiB；lil 有 1,775,504 个节点，逻辑节点约 27.09 MiB，容量同为 32 MiB。
结束导出时两端数组会短暂共存；公共结果键、Python 结果表和其他容器均另计。
因此这些索引大小不能与进程峰值 RSS 直接相减。

## 历史测量和范围

本轮修改前已按上轮 prefix 实验 JSON 的原 SHA-256 归档九个测量文件至
`validation/snapshots/prefix_trie_20261116_abi2`；更早 CASM 的 ABI-1 快照也保留。
旧原始数据不改写。当前默认 C++ 完整代表存储发生变化，旧脚本重新运行时须
明确选择存储模式，旧 CASM 时间不当作本轮实测。

本轮未重跑 CASM 外部程序，也未修改论文表 2、正式性能 CSV 或 CIF/DFT 管线。
论文算法表述可从统一说明抽取，性能表若切换新实现须按正式同群同输出口径重测。
