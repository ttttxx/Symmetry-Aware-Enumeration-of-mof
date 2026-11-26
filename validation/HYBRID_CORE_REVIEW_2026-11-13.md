# Python 子群调度 + C++ 搜索核心

## 本轮交付

作者接受完整字节键去重、保留子群 DFS，并提出混合 C++ 核心。本轮实现可构建、可回退、可显式选择的后端：

- `native/enumeration_core.cpp`：C++17 C ABI，DFS/计数器两种核、canonical 前缀比较、active 稳定子证书、多个 64 位字的覆盖掩码、完整字节键 `unordered_set`、数量上限及 steady-clock 搜索截止。
- `enumerator_native.py`：`HybridMOFEnumerator`，用标准库 ctypes 调用。Python 保留群/子群任务组织、排除条件编译、输出字典、原标签解码、排序与后续 CIF 构建。
- 每个子群任务调用一次核心，结果成批复制回 Python；没有逐搜索节点 Python 回调。原生上下文跨同次枚举的任务保留完整键去重，结束/异常时释放。
- `tools/build_enumeration_core.py`：macOS/Linux C++17 构建，无 pybind11/nanobind 或额外 pip 依赖。当前已在本机 Apple Clang 21 arm64 编译通过；Linux 构建分支提供但本轮未在 Linux 执行。共享库忽略 Git，目标机器需重建。
- `main.py --enumeration-backend {legacy,python,auto,cpp}`，默认仍为 legacy。python 强制优化字节键 Python DFS；auto 可用时使用 C++、不可用/超出支持维度时回退；cpp 强制要求共享库可用。Python `run_enumeration` API 增加末尾可选 `enumeration_backend` 参数，保留默认及现有位置参数。
- CLI benchmark JSON 加入 requested/active backend；策略名区分 byte-key subgroup DFS 和 legacy。与 exact-rho 配额组合会在参数解析阶段明确拒绝，尚未给配额子类复用整个 Fix(H) 覆盖证书。

构建/使用/边界说明见 `native/README.md`，主 README 已增加入口。

## 正确性与边界

- 同一标签字典序编码、原池顺序、子群任务顺序、全群规范代表、稳定子、轨道大小、完整结果顺序和 ID 均与字节键 Python DFS 对照。
- 原生字节标签支持至多 256 种、1–512 位点；auto 对不支持范围走 Python，保留其完整 tuple 回退。不是通过截断键降低数据精度。
- 全局输出 cap 和恒等任务新增唯一代表预算被转换为原生唯一结果数量上限。只有完整任务才登记覆盖；结果传回后仍由原有 Python 调度检查上限/时间状态，采用保守耗尽判定。
- task deadline 在参数准备后折算为剩余时长，核心使用 steady clock 协作检查。结果批量复制和 Python 记录物化可能在搜索截止之后结束；这是软截止，不承诺整个调用在某个精确毫秒返回。速度变化也不保证固定秒数返回相同集合。
- 编译核心保留完整键碰撞检查。Python 侧也保留成员索引，支持既有查重接口；C++ 内存和 Python 输出在枚举阶段同时存在，本轮没有报告新的独立进程内存基准。
- 无导入时编译或安装操作。显式 cpp 失败会给出构建提示；运行中内核错误会抛出，不中途静默重跑任务。

## 验证

```
.venv/bin/python tools/build_enumeration_core.py
clang++ -std=c++17 -Wall -Wextra -Werror -fsyntax-only native/enumeration_core.cpp
.venv/bin/python -m pytest -q tests/test_native_enumerator.py tests/test_enumerator_storage.py tests/test_counter_enumerator.py tests/test_subgroup_task_dedup.py tests/test_identity_budget.py
.venv/bin/python main.py --help
.venv/bin/python validation/benchmark_hybrid_20261113.py
```

- 本机共享库构建、Clang 严格语法/警告检查和 CLI help 均通过。
- 64 tests passed，19 条既有 spglib 弃用提示。新增覆盖 C++ DFS/计数器与 Python 完整输出、总量/新增预算、缺库 auto 回退、强制 cpp 报错、257 标签回退、输入池顺序、重用上下文、零上限、tuple 索引选项、超过 64 位的排除掩码、截止前/中途超时不登记整域覆盖，以及 CLI/default 选择和 rho 互斥。
- 五个真实算例各三版本三次穷尽执行，45 次结果列表逐项完全相同，并检查整数轨道大小和等于 |X|；另外五例 × 三种数量设置（总上限 7/100、恒等新增配额 3），15 组均比较三版本完整返回列表相同。
- 新核心计数器和 DFS 使用相同 C++ 共享库、相同子群调度、同一 exact-key 成员结构、同样输出物化。满足“同语言/同群/同输出工作量下再比较”的约定。
- 未重跑完整 CIF/GCMC/DFT 管线；未改论文及正式比较 CSV。默认运行仍走旧入口，新增后端需要显式选择。

## 同轮热运行结果

同一进程，预导入 SymPy，各版本三次，轮换顺序，每次计时前收集垃圾，计时内保持默认 GC。时间包含 Python 子群构建、覆盖编译、接口编解码、原生/解释器搜索、结果物化、排序和验证，不包含图/群分析及共享库构建。

| 算例 C3/C4 | Python 字节 DFS / s | 混合 C++ DFS / s | 混合 C++ 计数器 / s |
|---|---:|---:|---:|
| ffc 2/2 | 0.001578 | 0.001780 | 0.001770 |
| ffc 3/3 | 0.002673 | 0.002176 | 0.002290 |
| fjh 1/2 | 0.041974 | 0.034876 | 0.034699 |
| fjh 2/1 | 0.052333 | 0.041593 | 0.048103 |
| fjh 1/3 | 0.366369 | 0.095167 | 0.140422 |

完整原始数据：`HYBRID_BENCHMARK_2026-11-13.json`。最大例端到端约快 3.85 倍；最小例跨语言固定开销抵消收益，因此不能声称所有规模都加速。

最大例最后一轮核心调用累计约 DFS 0.0226 s / 计数器 0.0633 s，计时不含 Python 数组准备和结果物化；它是单轮诊断值，不当作三轮中位数。端到端混合版本计数器约比 DFS 慢 48%，支持继续子群 DFS。fjh 1/2 的两种混合端到端时间非常接近，受共同约 30 ms 子群准备支配，不解释为稳定算法反转。

## 下一步

1. 中大型枚举可通过显式参数使用混合 DFS；缓存 Python 子群/partition 构建是小例的下一优先项，首次和缓存命中计时分别报告。
2. 最大例原生核约 23 ms，剩余端到端成本包含约 30 ms 任务准备、约束编译、跨语言转换与 Python 结果物化。进一步优化先拆开这些阶段，不将全部余项说成 C++ 搜索慢。
3. 正式重测表 2 时统一群、输出工作量和冷/热边界，使用已记录的 active backend；不直接把本轮热运行结果覆盖旧论文数值。
