# main_MATCH 当前算法重测（2026-11-16）

本轮只生成独立结果与论文修改建议，不改论文和旧 CASM 数据。

```sh
.venv/bin/python experiment/match_refresh_compare.py
python3.11 experiment/match_vba_refresh.py
```

穷尽对照每侧三次独立顺序进程；包括当前 cpp/trie/prefix、CASM 同 Γ 同结果、同 Γ 只计数、原生 F 只计数，另外五例运行原 legacy。共 123 次。

当前后端强制 ABI 3、native_store=trie、prefix_trie=True；task_time_limit=None、max_total_results=None、identity_orbit_limit=inf。双方同 Python 3.11，图/群构建、导入和验证在计时外；子群任务准备、保存和排序在计时内。CASM 指定实际置换集合一致的 Γ。验证沿用独立 Burnside、整数覆盖、归一完整轨道记录摘要。

穷尽 harness 沿用 JSONL 断点续接。要重新测量请指定全新路径，例如：

```sh
.venv/bin/python experiment/match_refresh_compare.py --output validation/MATCH_REFRESH_new.json
```

vba 单列普通子群 DFS 旧/新对照，两类节点输出 41984、三类输出 60549，采用论文中的实际节点池，关闭超时，比较完整有序输出及稳定子分布。它不测精确 ρ 配额后端，也不进行 (2^{32}) 或 (3^{32}) 全空间 CASM 枚举。

结果与论文修改清单见 `validation/MATCH_ALGORITHM_PAPER_REVIEW_2026-11-16.md`。输入、源码、原生库哈希见 `validation/MATCH_REFRESH_MANIFEST_2026-11-16.json`，冻结版本见其 snapshot_root。默认 worker 为 python3.11，CASM 2.3.1；生产默认后端未切换。
