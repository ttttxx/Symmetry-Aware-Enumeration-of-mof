"""Aggregate recorded measurements only; do not infer unmeasured counters."""
import json
from pathlib import Path
import statistics
ROOT=Path(__file__).resolve().parents[1]
V=ROOT/'validation'
j=json.loads((V/'MATCH_REFRESH_2026-11-16.json').read_text())
v=json.loads((V/'MATCH_VBA_REFRESH_2026-11-16.json').read_text())
rows=[]
for s in j['summary']:
    rr=[r for r in j['records'] if (r['topology'],r['c3'],r['c4'])==(s['topology'],s['c3'],s['c4']) and r['side']=='ours']
    assert len(rr)==3
    stats=rr[0]['task_stats'];index=rr[0]['index_stats']
    for r in rr:
        for k in ['tasks','completed_tasks','search_nodes','leaf_visits','duplicate_leaves','pruned_prefixes']:
            assert r['task_stats'][k]==stats[k]
        assert r['index_stats']==index
    rows.append(dict(case=f"{s['topology']} {s['c3']}/{s['c4']}",mode='exhaustive',sites=rr[0]['sites'],group=rr[0]['group_size'],naive=s['naive'],returned=s['orbits'],
        **{k:stats[k] for k in ['tasks','completed_tasks','search_nodes','leaf_visits','duplicate_leaves']},
        completed_region_prunes=stats['pruned_prefixes'],canonical_rejections=None,label_choice_attempts=None,old_dfs_nodes=None,casm_dfs_nodes=None,
        prefix_hits=rr[0]['prefix_stats']['hits'],casm_candidates=s['naive'],
        legacy_s=s.get('legacy_median_s'),ours_s=s['ours_median_s'],casm_records_s=s['casm_records_median_s'],casm_count_s=s['casm_count_median_s'],
        kernel_s=statistics.median(r['native_kernel_s'] for r in rr),coverage_compile_s=statistics.median(r['task_stats']['compile_s'] for r in rr),
        legacy_rss_mib=s.get('legacy_median_rss_mib'),ours_rss_mib=s['ours_median_rss_mib'],casm_rss_mib=s['casm_records_median_rss_mib'],
        legacy_speedup=s['legacy_median_s']/s['ours_median_s'] if 'legacy_median_s' in s else None,casm_speedup=s['matched_records_speedup'],
        trie_nodes=index['nodes'],trie_logical_mib=index['node_bytes']/1048576,trie_capacity_mib=index['capacity_bytes']/1048576,record_checks='Burnside + integer coverage + normalized complete records'))
for s in v['summary']:
    rr=[r for r in v['records'] if r['c4']==s['c4'] and r['side']=='ours']
    stats=rr[0]['task_stats']
    for r in rr:
        for k in ['tasks','completed_tasks','search_nodes','leaf_visits','duplicate_leaves','pruned_prefixes']:
            assert r['task_stats'][k]==stats[k]
    rows.append(dict(case=f"vba C4={s['c4']}",mode='bounded',sites=32,group=rr[0]['group_size'],naive=s['naive'],returned=s['returned'],
        **{k:stats[k] for k in ['tasks','completed_tasks','search_nodes','leaf_visits','duplicate_leaves']},
        completed_region_prunes=stats['pruned_prefixes'],canonical_rejections=None,label_choice_attempts=None,old_dfs_nodes=None,casm_dfs_nodes=None,
        prefix_hits=rr[0]['prefix_stats']['hits'],casm_candidates=None,
        legacy_s=s['legacy_median_s'],ours_s=s['ours_median_s'],casm_records_s=None,casm_count_s=None,
        kernel_s=statistics.median(r['kernel_s'] for r in rr),coverage_compile_s=statistics.median(r['task_stats']['compile_s'] for r in rr),
        legacy_rss_mib=s['legacy_median_rss_mib'],ours_rss_mib=s['ours_median_rss_mib'],casm_rss_mib=None,
        legacy_speedup=s['legacy_speedup'],casm_speedup=None,trie_nodes=None,trie_logical_mib=None,trie_capacity_mib=None,
        record_checks='all ordered records/IDs agree; bounded, not full-space coverage'))
def num(x):return '—' if x is None else (f'{x:,}' if isinstance(x,int) else f'{x:.6g}')
def pair(*xs):return ' / '.join(num(x) for x in xs)
headers=['算例','模式','位点/群阶','原始空间 N','输出 r','完成/调用任务','DFS 状态','完整叶/CASM 完整候选','覆盖剪枝','重复叶/前缀命中','时间：旧/新/CASM（s）','速度比：旧/新、CASM/新','RSS：旧/新/CASM（MiB）','C++核（s）','覆盖编译（s）','CASM只计数（s）','结果Trie节点','Trie逻辑/容量（MiB）']
md='''# 枚举算法多指标对比总表（2026-11-16）

下表仅汇总本轮已经记录的数据，不重新计时。时间和 RSS 均为三次独立进程的中位数；计数在三次运行中一致。CASM 主列指相同 Γ、保存相同结果字段的配置。`—` 表示未测或不适用，不能解释为零。旧版指原 Python legacy；新版指 ABI-3 C++ 子群 DFS、完整代表压缩 Trie，显式开启前缀覆盖。

'''
md+='| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'
for r in rows:
    cells=[r['case'],'穷尽' if r['mode']=='exhaustive' else '限量',pair(r['sites'],r['group']),num(r['naive']),num(r['returned']),pair(r['completed_tasks'],r['tasks']),num(r['search_nodes']),pair(r['leaf_visits'],r['casm_candidates']),num(r['completed_region_prunes']),pair(r['duplicate_leaves'],r['prefix_hits']),pair(r['legacy_s'],r['ours_s'],r['casm_records_s']),pair(r['legacy_speedup'],r['casm_speedup']),pair(r['legacy_rss_mib'],r['ours_rss_mib'],r['casm_rss_mib']),num(r['kernel_s']),num(r['coverage_compile_s']),num(r['casm_count_s']),num(r['trie_nodes']),pair(r['trie_logical_mib'],r['trie_capacity_mib'])]
    md+='| '+' | '.join(cells)+' |\n'
md+='''
## 指标含义与不能混用的地方

- **DFS 状态**：C++ `Task::dfs` 的进入次数，包含中间状态、叶层、已覆盖区域入口及中断检查，跨任务累计；在规范比较处被拒绝、未进入递归的标签尝试不计入。因此它既不等于完整候选数，也不等于全部比较次数。CASM 本轮路径没有 DFS，CASM 的 DFS 状态数不适用。
- **完整叶/CASM 完整候选**：前者是本文到达 `accept` 的完整规范叶次数，后者是 CASM 计数器生成并检查的全部完整指派。二者位于不同算法阶段，不能直接据其比例声称总操作量减少。所有本轮本文完整叶均为新代表，重复叶为 0。
- **覆盖剪枝**：原生 `Task::excluded` 中已完成任务区域等式证书命中的次数。它对应原始 JSON 的 `pruned_prefixes`，**不包含全部规范前缀拒绝数**，不能拿它除以节点数报告总剪枝率。
- **前缀命中**：已完成搜索前缀缓存命中数，全部为 0。vba 在最后的总上限处保存证书，但没有后续任务继续利用；这部分不能计为实测加速来源。
- **C++核**：跨任务累计的原生调用时间，包含该任务内搜索、Trie操作等。覆盖编译是单列子阶段，不等于全部子群准备。各阶段分属独立计时区间；中位数之和不必恰好等于总耗时中位数，不能把总时间减去核时间的余项全部叫作 Python 或排序成本。
- **Trie逻辑/容量**：结果索引节点的逻辑空间与原生 vector 保留空间，不含结果键、Python结果容器、解释器等；原生索引导出期间还可能暂存双份。vba worker 未记录索引节点统计，不能用返回数推算填补。
- **RSS**：进程截至枚举结束的峰值，含依赖导入/输入准备产生的内存，不是单独搜索工作区，也不是整机内存占比。
- vba 两行限定输出 41,984 / 60,549；这些行不代表穷尽，更没有对应的 CASM 全空间计时。vba 的完整有序结果、ID及元数据与旧版相同；其余九行通过独立 Burnside、整数轨道覆盖和完整归一记录校验。

## 时间与空间复杂度

记 n 为位点数，g=|Γ|，q 为最大标签池大小，N=∏v|A_v|，r 为输出轨道数。以下分析对应保存完整记录的穷尽模式，公共图/群输入准备单列。

旧版 Python 与新版后端的穷尽搜索最坏情况均为指数级；此次优化没有改变这个结论。旧版未压缩完整键 Trie 最坏需要 O(nr) 个节点；当前压缩 Trie 为 O(r) 个节点，但双方完整输出仍占 O(nr) 空间，因此不能说总程序空间复杂度从 O(nr) 降为 O(r)。

**CASM 本轮完整候选路径**必须生成 N 个完整指派。一次规范性/稳定子检查最坏扫描 g 个操作和 n 个位点；含结果保存及排序的粗上界为 `O(N g n + n r + r log r)`，即标签数固定且 n 增长时仍为指数级。短路比较可以降低实际工作量；该上界不是其每次运行的实测成本。

**本文子群 DFS**仍为最坏指数级。对任务 H 的 m_H 个变量，设按顺序的域大小为 q_H,j；忽略剪枝时，该任务搜索树的前缀数量为 `S_H = ∑(k=0…m_H) ∏(j=1…k) q_H,j`，空乘积为 1。多任务的结构性搜索上界随 `∑H S_H` 增长；恒等任务本身包含完整空间 N。任务准备、群比较、覆盖检查和结果存储仍要额外计费，不能把运行时间直接写成 `O(r)`，也不能把子群格或 Trie 解释为多项式复杂度保证。

更具体地，令 A 为全部标签选择尝试数（包括规范比较后拒绝而未进入 DFS 的尝试），L 为进入完整叶处理次数，C_cov 为覆盖检查总工作量，T_task 为任务准备成本，则当前实现的一个保守成本表达为：

`T_ours = O(T_task + A g n + C_cov + L n q + n r + r log r)`。

其中 `O(nq)` 是当前排序兄弟链紧凑 Trie 的最坏单次插入粗界；q≤256，固定标签数时为 O(n)。DFS 游标复用减少实际重复比较，但本轮没有记录 A、群比较次数或全部规范拒绝数，**不能用已记录的 DFS 状态数替代 A 代入复杂度公式**。该表达分开列出了工作来源，不是给出紧确渐近界。

双方保存 r 个长度 n 的完整代表，都至少需要 Ω(nr) 的输出时间与存储。当前结果键/结果表空间为 O(nr)，本文压缩结果 Trie 另有 O(r) 节点，逻辑节点数至多 2r+1；实际总空间还包括输入群、子群划分、DFS状态、覆盖条件、可选前缀缓存与跨语言导出缓冲。Trie节点的 O(r) 不能写成整个程序的 O(r) 内存。若讨论只计数模式，应另列结果保留需求，不能直接套用本轮完整结果 RSS。

## 尚未记录、不能从现表推出的指标

旧 legacy 的 DFS状态数、全部规范前缀拒绝数、标签尝试次数、实际群比较次数、各类剪枝耗时、任务构造与结果解码/封送/排序的完整分项、只计搜索阶段的增量峰值内存，本轮均未完整采集。CASM 已验证完整候选访问数为 N，但其内部逐操作比较数未记录。若论文要解释具体每项改进贡献，需另做插桩/消融；插桩计时应与正式性能计时分开。

汇总脚本：`validation/summarize_match_metrics_20261116.py`。输入为 `MATCH_REFRESH_2026-11-16.json` 与 `MATCH_VBA_REFRESH_2026-11-16.json`。结构化总表为 `MATCH_METRICS_COMPARISON_2026-11-16.json`。本次未修改论文、算法或原始计时文件。
'''
(V/'MATCH_METRICS_COMPARISON_2026-11-16.md').write_text(md)
(V/'MATCH_METRICS_COMPARISON_2026-11-16.json').write_text(json.dumps(dict(rows=rows,note='Measured fields only; missing counters are null; theory explained in accompanying Markdown.'),indent=2)+'\n')
print(f'Wrote unified table: {len(rows)} cases, {len(headers)} columns; counters consistent across repetitions.')
