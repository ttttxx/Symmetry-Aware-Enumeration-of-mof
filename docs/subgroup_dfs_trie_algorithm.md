# 子群固定点枚举：规范前缀 DFS 与 Trie

对应实现：`enumerator.py` 的任务调度、`enumerator_task_dedup.py` 的任务覆盖编译、
`enumerator_native.py` 的任务接口和 `native/enumeration_core.cpp` 的 C++ 搜索核。
本文描述当前 C++ Trie 方案，覆盖前缀复用开启时的完整流程。CLI：

```sh
.venv/bin/python tools/build_enumeration_core.py
.venv/bin/python main.py --enumeration-backend cpp --native-store trie --prefix-trie
```

## 1. 输入、输出与规范代表

输入为顶点集 $V=\{0,\ldots,n-1\}$、各顶点允许标签集 $A_v$，以及作用于顶点的
有限置换群 $\Gamma$。要求群作用保持允许标签集；当前项目通过节点度数规定
标签池，所用置换保持节点度数。标签按统一顺序编码为整数。

候选空间和群作用为

$$
X=\prod_{v\in V} A_v,\qquad (g\cdot x)_v=x_{g^{-1}(v)}.
$$

每个轨道选择字典序最小的指派

$$
\operatorname{can}(x)=\min_{g\in\Gamma}g\cdot x.
$$

输出记录包含规范指派 $x$、稳定子大小 $s(x)=|\{g:g\cdot x=x\}|$，以及
轨道大小 $|\Gamma|/s(x)$。最终按稳定子大小降序稳定排序，保留高对称优先的
结果顺序。返回记录的 `id` 在排序后重编号；内部 `classes` 键和 Trie 的样本
记录下标始终使用插入顺序，二者须区分。

关闭预算并使全部任务搜索结束时，每个轨道返回一个代表。受限运行返回该顺序下
已经发现的部分代表。本实现枚举给定标签池的全部相容指派，未增加固定组成计数
约束，精确稳定子层配额现通过独立调度接入相同优化搜索核，见第 9 节。

## 2. 子群任务及搜索变量

子群 $H\leq\Gamma$ 将顶点划分为轨道
$P_H=(B_0,\ldots,B_{m-1})$。任务搜索的固定点区域为

$$
\operatorname{Fix}(H)=\{x\in X:x_u=x_v\text{ whenever }u,v\in B_i\}.
$$

DFS 每次给整个 $B_i$ 赋同一个标签，形成 $m$ 个任务变量。允许池在每个块内
相同，当前代码由块的首个顶点取得该池。初始任务包含全群；其余子群由 SymPy
构造，合并产生完全相同顶点划分的任务，并加入恒等任务。任务按变量数从少到多
执行，同分保持生成顺序。恒等任务的每个块只含一个顶点，其搜索域为整个 $X$。

不同子群的固定点区域可以相交。算法在搜索时复用已完成区域，并由完整代表 Trie
保证全局唯一。子群优先次序组织搜索、帮助较早发现高对称指派；完整覆盖由恒等
任务及安全剪枝共同保证。

## 3. 搜索状态和两种 Trie 记录

一次枚举维护结果表 $R$、完整代表 Trie $T_R$、完整任务覆盖表 $C$，以及可选的
前缀覆盖 Trie 集合 $T_C$。各任务共享这些状态。

| 状态 | 记录内容 | 可以做什么 |
|---|---|---|
| $R$ | 已输出规范指派及稳定子/轨道大小 | 返回和构建模型 |
| $T_R$ | 长度为 $n$ 的完整规范键 | 精确查重；首次出现时加入 $R$ |
| $C$ | 已完整结束任务的顶点划分 | 排除整个已完成固定点区域 |
| $T_C$ | 来源任务划分、允许标签域、已完成的变量前缀 | 约束包含关系成立时跳过子树 |

两种 Trie 记录的事实不同。$T_R$ 的终点表示一个完整指派已经保存；$T_C$ 的
关闭节点表示来源任务在这一前缀下的规范搜索已全部结束。完整代表出现过一次，
不能据此关闭其所在的搜索前缀。逻辑上分别维护这两类标记，避免把查重当作覆盖。

DFS 状态包含当前部分指派、当前变量深度、仍可能与指派相等的群操作及其比较
游标、整域覆盖条件的活动掩码，以及各层已完成选项数。标签采用全局整数顺序，
任务变量次序与结果键的顶点次序分别保留。

## 4. 一个任务中的 DFS

### 群作用前缀比较

对每个非恒等操作 $g$，维护已经确认相等的顶点前缀长度 $c_g$，比较
$x_v$ 与 $x_{g^{-1}(v)}$。从 $c_g$ 向后推进：

1. 任一端尚未赋值时，保留该操作和游标，等待更多变量确定。
2. 首次差异满足 $x_{g^{-1}(v)}<x_v$ 时，整棵当前子树可剪。
3. 首次差异满足 $x_{g^{-1}(v)}>x_v$ 时，该操作的规范条件已满足，后代无需再比。
4. 比到第 $n$ 位仍相等时，该操作在这一完整指派的稳定子中。

游标只越过已确定且相等的位置，首次确定差异在后续赋值中保持不变。因此不会
因尚未赋值的位置误剪。到完整叶子时，所有保留操作都与当前指派相等，稳定子
大小直接为 $1+|\text{active}|$，其中 1 为恒等操作。

### 完整任务区域排除

对于 $C$ 中的来源划分 $P$，将块内等色条件编译到当前任务变量上。相同变量
无需检查，单标签域可作为常量，其余条件在对应变量确定时执行。活动位保留尚未
被违反的来源区域；若某一活动区域的全部等式已经确定，当前所有补全均落入该
已完成区域，可以立即返回“该子树已处理”。掩码使用多个 64 位字。

### 中断任务的前缀区域排除

任务中断时，DFS 各层已完成选项数与中断路径一起确定已完成的兄弟子树。只把
这些子树写入来源任务的 $T_C$，中断路径保持开放。完整完成的任务直接加入
$C$，无需逐前缀建树。

后续任务的某个分支只有同时满足下列条件，才可使用来源 Trie：

- 当前任务各顶点允许标签集均包含于来源任务的允许集。
- 当前分支的**所有补全**都满足来源任务的块内等色约束。未赋值顶点若属于同一
  当前变量，可证明其相等；不同未赋值变量不能默认相等。
- 按来源任务的变量次序查询，已确定标签到达一个关闭前缀。

代码预计算这些条件和前缀何时首次可判定，只在相应 DFS 深度查询。标签域包含
检查采取保守判据，允许少剪，不会把尚未覆盖区域剪掉。缓存容量耗尽时停止增加
相应节点，继续正常搜索。

### 完整叶子的处理

通过上述检查的完整叶子已经是规范代表，无需再构造完整群轨道或重复规范化。
在 $T_R$ 中执行精确插入：若完整键已有终点则跳过；若首次插入则追加结果与
由 active 集合得到的稳定子。重复键不消耗新增结果预算。

## 5. 统一伪代码

`CLOSED` 表示当前任务约束下该子树的规范搜索已经处理完成，可以复用；`OPEN`
表示预算或时间中断。被安全剪枝的子树也可视为处理完成，因为其中没有尚待输出
的规范叶子。预算检查与当前实现一样放在进入 DFS 及每轮选项之前。

```text
ENUMERATE(Γ, {A_v}, budgets):
    encode labels in a common lexical order
    R ← empty result table
    T_R ← empty compact complete-key Trie
    C ← empty complete-task coverage table
    T_C ← empty completed-prefix coverage cache
    tasks ← subgroup fixed-point partitions, ordered by variable count

    for P in tasks:
        if global output limit reached: break
        compile C into equality events and region masks for P
        prepare eligible T_C queries and their decision depths
        initialize coloring, group cursors and DFS progress
        state ← DFS(P, depth=0)

        if state=CLOSED and task did not trigger a cutoff:
            add P to C
        else if state=OPEN and prefix caching is enabled:
            compress completed siblings from the saved frontier into T_C

    stable-sort R by decreasing stabilizer size; renumber public output IDs
    decode labels and return R

DFS(P, k):
    if a search budget is reached:
        save the first interrupted frontier
        return OPEN
    if a completed task region covers every completion: return CLOSED
    if a compatible closed prefix covers every completion: return CLOSED

    if k = number of blocks in P:
        s ← 1 + number of surviving nonidentity operations
        if T_R.insert(full coloring) reports a new key:
            append (coloring, s, |Γ|/s) to R
        return CLOSED

    progress[k] ← 0
    for label in allowed choices for block P[k]:
        if a search budget is reached:
            save the first interrupted frontier
            return OPEN
        assign label to all vertices in P[k]
        advance the group-comparison cursors
        child_state ← CLOSED if infeasible or canonically rejected
                      otherwise DFS(P, k+1)
        restore assignment and comparison state
        if child_state=OPEN: return OPEN
        progress[k] ← progress[k]+1
    return CLOSED
```

实现中 C++ 执行整个任务的 DFS；Python 负责划分准备、预算参数、每任务批量取
结果及最终接口。没有逐搜索节点的 Python 回调。C++ 完整代表 Trie 跨任务存活，
在上下文释放前导出节点数组，以保持枚举后的查询能力。

## 6. 完整代表 Trie 的紧凑数组实现

`CompactTrie` 是路径压缩的字节 Trie。节点连续存入 `vector<Node>`，每个节点
固定 16 字节，四个字段均为 32 位无符号整数：

```text
Node = (sample_record, depth, first_child, next_sibling)
```

父节点深度为 $a$、当前节点深度为 $b$ 时，边标签为
`R[sample_record].key[a:b]`。结果键在原生结果表中保存一份，Trie 不再复制边
字符串。数组扩容后整数下标仍有效。兄弟节点按首字节排序，0 和 128–255 都按
无符号字节比较，空字节不作为终止标记。

插入时沿边比较完整标签区间：

1. 缺少对应首字节分支时，增加一个指向完整键末端的叶子。
2. 在压缩边内部首次不同时，在最长公共前缀处拆边，增加分叉节点与新叶子。
3. 匹配到长度 $n$ 的已有叶子时，判为重复。

所有键长度相同，叶子即终点。每个非根内部节点至少有两个孩子，存储 $r$ 个
不同键所需节点数至多 $2r+1$。节点逻辑空间为 $16N$ 字节；`vector` 实际保留
容量另记，结果键和公共输出表的空间另计。当前用 32 位下标，超过索引范围会
明确报错，不能静默截断。

Trie 模式不建立 C++ 完整键 `unordered_set`，也不建立 Python 完整键哈希副本。
Python 保留按整数下标访问的结果表，以及一份导出的连续节点数组；查询时仍沿
Trie 精确比较，边标签从相应结果记录读取。原生上下文销毁后不保留悬空指针。

目前仅完整代表 Trie 使用上述压缩数组；前缀覆盖缓存保留带来源任务约束的独立
树结构，其记录上限和生命周期独立。这不影响统一 DFS 中各项检查的先后关系。

## 7. 正确性与成本边界

规范前缀剪枝只丢弃已证明字典序更大的补全。整域覆盖和关闭前缀只丢弃此前已
完整处理、或已由其他安全条件排除的规范搜索区域。完整代表 Trie 精确比较整个
键，保证每个保存的规范指派只有一份。恒等任务覆盖 $X$，因此穷尽运行会访问或
安全排除每个规范指派。由有限群作用的唯一字典序最小元，输出与 $X/\Gamma$
一一对应。

完整键 Trie 的一次查询/插入比较至多 $n$ 个边标签字节，另有兄弟分支查找。
当前兄弟链最坏每层扫描至多 $q\leq256$ 个标签，粗界为 $O(nq)$；空间节点数
为 $O(r)$，公共结果键仍占 $O(nr)$。字节哈希的期望查询成本包含 $O(n)$ 的键
哈希与可能的键比较。二者实际速度受标签数、共享路径、分配成本和缓存行为影响。

总搜索仍可能指数增长；恒等任务最坏涉及 $\prod_v|A_v|$ 个指派。子群构造、群
作用比较、覆盖条件查询及输出物化均计入端到端成本。Trie 替换完整代表哈希会
改变索引成本，不会自动减少 DFS 节点数。相关本轮实测见
`validation/COMPACT_TRIE_REVIEW_2026-11-16.md`，不以存储结构的名字代替性能证据。

## 8. 实现与运行边界

| 算法部分 | 当前代码入口 |
|---|---|
| 子群/任务顺序、输出排序 | `MOFEnumerator.enumerate_subgroups` / `enumerate` |
| 整域等式事件、完成任务登记 | `TaskDedupMOFEnumerator` |
| 整数编码、Python 对照 DFS | `CasmInspiredMOFEnumerator` |
| 参数封送、结果批量导出、数组索引保留 | `HybridMOFEnumerator` |
| C++ 规范比较、覆盖查询、叶子接受 | `Task::dfs` / `prefix_covered` / `accept` |
| C++ 完整代表压缩数组 Trie | `CompactTrie::insert` / `lookup` |

当前 C++ 限制为 1–512 个顶点、至多 256 个编码标签，ABI 为 4。强制 `cpp`
加载失败时明确报错，`auto` 可回退 Python 字节/元组哈希 DFS，并报告实际后端。
主程序与API未选择后端时默认 auto，优先C++，不可用时回退优化Python；legacy须显式选择。`native_store="hash"` / `--native-store hash`
用于同核对照；C++ 默认采用 `trie`。前缀覆盖开关仍为显式可选项。

时间截止为协作截止，任务参数准备、前缀快照、结果导出及排序存在固定成本；
固定秒数运行不保证不同版本返回相同结果集合。精确稳定子层配额将覆盖证书限定在同一目标稳定子阶数内复用，见第 9 节。


## 9. 精确稳定子层配额（2026-11-23）

`HybridRhoQuotaMOFEnumerator` 复用原 `RhoQuotaScheduleMixin`：按子群阶数递减，
同阶合并相同划分并按原字典序执行任务，保留输入标签遍历顺序。Python 和 C++
均复用普通优化 DFS 的整数编码、比较游标和稳定子提取。完整叶的稳定子阶数
必须等于当前目标阶数 s 才允许插入全局结果索引；其他阶数不写入结果索引，
也不消耗当前层预算。单任务新增结果上限取该层剩余配额，并同时执行全局上限。

覆盖证书的语义是“该固定点域中目标稳定子阶数为 s 的规范搜索已完成”。因此，
整域覆盖表和中断前缀缓存只在同一目标层内共享；进入下一层时二者清空，
全局完整代表索引和已接受结果保留。只有没有触发超时/数量上限的完整任务
能登记该层整域证书。前缀缓存仍只保存已经完成的兄弟子树；开放路径保持开放。
配额填满仍报告 quota_reached，全部任务完成且不足配额才报告 layer_exhausted。

原 legacy 配额保留，用于回归与默认入口；其 deadline 参数改为关键词传递，
修复旧位置调用将截止时刻传入 identity_budget 的问题。新 CLI/API factory
支持 python/auto/cpp 与 rho_quota_per_level 同时指定，ABI 已升至 4，须重建库。
历史 ABI-3 普通枚举数据与源码快照保持原样。当前真实配额重测和验证见
`validation/RHO_NATIVE_VBA_2026-11-23.json` 与对应报告。
