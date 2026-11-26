# CASM 2.3.1 枚举实现与表 2 性能差距

## 范围与结论

本轮只做源码追查、性能诊断和消融；未修改 `enumerator.py`、正式 benchmark CSV 或论文。
对应主稿 `paper/main_MATCH.tex` 的 `tab:casm`。检查本机安装的 libcasm-configuration 2.3.1，并获取官方 v2.3.1 源码，commit `eec3bb720fa0d032f024cabe04bf11d7344db666`，临时检出位置 `/private/tmp/mof-casm-configuration-2.3.1`。该临时目录不是长期证据入口，下列固定 commit 链接可追溯源代码。

表 2 的旧时间比依次为约 661、60、59、5.2、5.3。小算例受首次导入和子群任务构建影响很大；最大算例的主要成本是 Python DFS 和叶节点规范化/稳定子/Trie 存储。减少子群重叠有价值，但不能解释或消除全部差距。

## CASM 实际调用链

`experiment/_casm_worker.py` → `ConfigEnumAllOccupations.by_supercell_list` → Python `_by_site` → C++ `ConfigEnumAllOccupationsBase` → primitive 检查 → C++ canonical 检查 → yield。

- [Python 枚举包装](https://github.com/prisms-center/CASMcode_configuration/blob/eec3bb720fa0d032f024cabe04bf11d7344db666/python/libcasm/enumerate/_ConfigEnumAllOccupations.py)：`_by_site` 的 while 循环逐候选调用筛选；Python 也在每个原始候选上运行，不只是最终输出。
- [C++ 候选计数器](https://github.com/prisms-center/CASMcode_configuration/blob/eec3bb720fa0d032f024cabe04bf11d7344db666/src/casm/configuration/enumeration/ConfigEnumAllOccupations.cc#L35)：`Counter<std::vector<int>>` 枚举每个位点允许的整数占位，`advance()` 更新复用的 `m_current`。每次仍会写入所枚举的全部位点，不能描述成只更新进位位点。
- [规范性判断](https://github.com/prisms-center/CASMcode_configuration/blob/eec3bb720fa0d032f024cabe04bf11d7344db666/include/casm/configuration/canonical_form.hh)：`is_canonical` 使用 `std::none_of`，遇到一个使候选更大的群操作就返回失败。CASM 取字典序最大，我们取最小；两种约定都可用于判等。
- [占位比较](https://github.com/prisms-center/CASMcode_configuration/blob/eec3bb720fa0d032f024cabe04bf11d7344db666/include/casm/configuration/ConfigDoFIsEquivalent.hh#L20)：各向同性占位以 `Eigen::VectorXi` 保存，通过 `A.permute_index(i)` 直接读取变换后位点值，首个不等位点决定大小。不为每个群操作分配完整变换构型。
- [置换索引](https://github.com/prisms-center/CASMcode_configuration/blob/eec3bb720fa0d032f024cabe04bf11d7344db666/src/casm/configuration/SupercellSymOp.cc#L130)：复用 supercell 的因子群/平移置换信息。
- [Python 绑定](https://github.com/prisms-center/CASMcode_configuration/blob/eec3bb720fa0d032f024cabe04bf11d7344db666/python/src/enumerate.cpp#L390)：`value()` 使用 `reference_internal` 返回当前构型引用。调用方若保留不同构型，应显式 `.copy()`；原 worker 仅计数，没有保存所有构型。

**这条具体路径没有子群格构建、子群固定点枚举或 DFS 前缀剪枝。** 本轮诊断直接计数原始访问数，五例分别为 128、2187、4096、65536、531441，与完整指派空间一致。不能把 CASM 其他模块/历史版本的算法推广到表 2 所用接口。

## 我们的额外工作

1. `enumerator.py:258` 在计时区间内首次导入 SymPy、构造子群并反复计算闭包/阶数，再将相同顶点轨道划分合并。它消除的是相同 partition，仍保留不同 partition 的固定点集合交集。全群任务也可能与随后子群任务形成相同 partition。
2. 最终恒等子群任务已经覆盖全部指派；穷尽模式下先前任务在覆盖意义上可以省去。受限输出时，这些任务承担优先获取高对称代表的作用，删除会改变输出策略。
3. `enumerator.py:162` 每层为 active 操作从位点 0 重扫已知前缀，创建下一层 active 列表；单选位点仍递归并参与比较。`fjh` 的 28 个位点中，两个二元/三元池算例分别只有部分位点可变。
4. `enumerator.py:144` 的叶节点已经过前缀比较，仍调用 `_get_canonical` 遍历全部群操作、每次构造完整字符串元组。后续可研究复用 DFS 的 canonical 不变量与 active 集合；须先验证该不变量和稳定子计数，不能直接全局删掉 `_get_canonical`。
5. `enumerator.py:97` 又遍历全群构造 `orbit_set` 并转 list；普通枚举器及分层配额子类的调用方均立即丢弃轨道列表，只使用稳定子大小。这里有直接可移除的中间分配。
6. Trie 只用于完整 canonical 键查询和插入，没有为 DFS 做前缀命中剪枝。最大算例 cProfile 记录约 646,700 次 `_TrieNode` 构造；每个节点还有 children 字典。固定长度键可优先比较 tuple/bytes 哈希表。
7. 我们保存 coloring、稳定子、轨道大小和全部记录，最后排序、检查覆盖；CASM 原 worker 只累计数量。两边输出工作量不同，但本轮 CASM 保存独立副本的诊断未显示显著时间增加，不能把全部差距归因于保存结果。

## 同机消融：保留现有 DFS 和叶节点实现

命令：`.venv/bin/python validation/profile_casm_gap_20261113.py`。
原始记录：`CASM_GAP_PROFILE_2026-11-13.json`。
每例独立进程；先执行首次 baseline，再取 3 次热运行中位数；identity-only 仅替换任务组织，保留现有前缀剪枝、canonical、Trie、稳定子、排序和覆盖检查。下表是热运行时间，不能直接覆盖论文的独立进程旧测量。

| 算例 C3/C4 | 完整任务 / s | 仅恒等任务 / s | 加速 | 叶节点总数 → 唯一代表 | 任务数 |
|---|---:|---:|---:|---:|---:|
| ffc 2/2 | 0.002280 | 0.000417 | 5.5× | 154 → 42 | 8 |
| ffc 3/3 | 0.008664 | 0.003676 | 2.4× | 1296 → 513 | 8 |
| fjh 1/2 | 0.098413 | 0.020050 | 4.9× | 1218 → 394 | 33 |
| fjh 2/1 | 0.191621 | 0.120166 | 1.6× | 6642 → 4756 | 33 |
| fjh 1/3 | 1.879108 | 1.460726 | 1.3× | 46467 → 36369 | 33 |

五例逐项比较了 canonical coloring→(stabilizer_size, orbit_size) 的完整映射，相同；同时整数轨道大小之和等于 |X|。这支持穷尽模式的消融结论，不证明受限输出相同。

插桩热运行的首次 DFS 前准备耗时：ffc 约 1.1–1.2 ms，fjh 约 29 ms。独立进程首次运行耗时还包含 SymPy 导入，首次 ffc 2/2 在本轮为 0.931 s，其余小例约 0.110/0.193/0.289 s；首例明显受瞬时冷启动影响，应保留原始记录，正式重测需更多重复并区分进程冷启动、首次导入和任务构建。

最大 fjh 的插桩统计：规范化约 0.434 s、稳定子/轨道集合约 0.379 s、Trie 插入约 0.342 s、Trie 查询约 0.038 s。cProfile 的热运行内，全群 `_apply_perm` 共 1,325,376 次；叶节点处理累计约占 profile 中枚举总耗时的 54%。这些是带插桩/采样工具开销的定位数据，不能与未插桩耗时精确相减；冷启动、任务构建和重复搜索也不能混称一项开销。

## 同群复核与比较口径

命令：`python3.11 validation/inspect_casm_path_20261113.py`。
原始记录：`CASM_PATH_CHECK_2026-11-13.json`。

原生 `by_supercell_list` 无直接自定义群参数，但底层 `is_canonical_configuration(..., subgroup=...)` 支持自定义操作列表。本轮筛出 proper 操作，将其位点置换与我们群的逆置换集合逐项核对，五例完全相同；不只核对群阶数。

同群直接枚举得到 **42、513、394、4756、36369**，均与我们一致。旧 harness 从 ours 和 full-group count 反推出 A/C 后重建 proper count 是代数恒等式，不能单独构成 CASM 的独立验证；本轮直接枚举补上了该证据。

| 算例 C3/C4 | CASM 原生 F 仅计数 / s | CASM 原生 F 保留副本 / s | CASM 底层 proper 子群 / s |
|---|---:|---:|---:|
| ffc 2/2 | 0.000129 | 0.000125 | 0.000125 |
| ffc 3/3 | 0.001696 | 0.001835 | 0.002057 |
| fjh 1/2 | 0.003338 | 0.003446 | 0.005132 |
| fjh 2/1 | 0.059907 | 0.059634 | 0.085539 |
| fjh 1/3 | 0.410168 | 0.408001 | 0.619452 |

均为同进程三次中位数。proper 诊断还改变了 Python 循环组织、每次传入 subgroup，并省去 volume-1 下无意义的 primitive 检查；因此 0.619 s 不是隔离群选择单一因素的严格消融。`F` 原生路径与同群诊断都仍明显快于当前大算例实现，足以支持优化核心代码的必要性。CASM 最大算例保留副本前后约 0.41 s 的小差异落在噪声内，不能解读为复制加速。

## 建议优化顺序

1. **先减叶节点重复工作**：稳定子仅计数、不要生成被丢弃的轨道集合；以定长哈希键替换只做完整键查询的 Trie；验证后复用 DFS 已知的规范性和稳定子候选。
2. **压低热循环成本**：标签内部整数化、预存逆置换、复用比较前缀进度、压缩只有一种选择的位点；经过 Python 层消融后，再决定 C++/编译核心的边界。
3. **分清运行目标**：穷尽计数/完整输出考虑单次 canonical DFS；高对称优先/分层受限输出继续保留子群组织。子群任务按有效位点/允许标签约束后的固定点空间去重，缓存群与 partition 构建。
4. **研究子群间重复**：H⊆K 时 Fix(K)⊆Fix(H)，要搜索互斥的“精确稳定子层”，须处理更大稳定子和共轭关系；仅用排除 Trie 或跳过一个已见叶节点不能省掉此前全部子树。先实施已有证据支持的低风险优化，再评估精确层搜索复杂度。
5. **正式重测**：同时报告准备/搜索/结果处理阶段，统一作用群、保存工作与冷/热口径，保留完整端到端时间。不能仅把准备移出计时就声称算法加速；先不改表 2 数值和论文措辞。

向导师解释时可说：现有实现为高对称优先输出承担了子群任务构建与重叠搜索成本，小算例还受首次导入影响；大算例另有可量化的规范化、稳定子和存储分配成本。CASM 此接口靠整数表示、复用状态和 C++ 短路规范性检查实现更低的候选处理成本。
