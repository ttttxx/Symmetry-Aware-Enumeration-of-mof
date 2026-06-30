# MOF Enumeration Project State

Last updated: 2026-06-30

## Project Scope

项目根目录是 `~/PyCharmMiscProject/mof`。

论文位于 `paper/main.tex`，DFT 验证包位于 `dft_validation/`。

本项目用群论和图自同构做 MOF 构型枚举。周期拓扑被抽象为 quotient graph，构型是节点到 building block 的 coloring。对称群 `Gamma` 作用在节点上，构型的 Symmetry Index 定义为：

```text
rho(c) = |Stab(c)| / |Gamma|
```

枚举器按对称性从高到低输出 symmetry-distinct colorings，并用 Burnside identity 验证完整性。物理结构由 `assembler.py` 用 Kabsch alignment 把刚性 SBU/linker 放回晶体，得到 CIF。

## Core Code

- `main.py`: `load_graph`、topology enrichment、degree requirements。
- `symmetry.py`: `SymmetryEngine`，计算图自同构和对称群。
- `enumerator.py`: `MOFEnumerator`，枚举 symmetry-unique colorings。
- `assembler.py`: `StructureAssembler`，Kabsch-based physical reconstruction。
- `dft_validation/select_structures.py`: 为 DFT 选择结构。当前采用固定组分内比较 `rho` 的设计，避免 composition/local chemistry confound。
- `dft_validation/run_dft_validation.py`: 阶段式流程 `select -> sanitize -> DFT -> analyze`；已修复 `fjh` topology 不再硬编码成 `ffc`；Phase 4 支持按 topology 分批并发运行多个 CP2K single-point jobs。
- `dft_validation/cp2k_runner.py`: CP2K input/run/parse。当前用 `.ssmp` OpenMP 单文件 CP2K，传统对角化 + Fermi smearing + Broyden mixing，并有 per-job OOM/time guard。
- `dft_validation/analyze.py`: 先按 topology，再按 composition 计算 Spearman `rho_s`，输出 CSV，并分别输出 `rho_energy_ffc.png` 与 `rho_energy_fjh.png`。

## DFT Validation Design

旧 `ffc` 只有 7 节点，只产生 3 个离散 `rho` 值，且组分内样本太少。当前结论是不支持“对称性越高越稳定”的强断言，根因是体系太小，能量主要被局部化学支配。

新方案引入 `fjh`：

- 28 节点。
- `|Gamma| = 16`。
- 枚举 394 个结构。
- `rho` 有 5 档：`1.0 / 0.5 / 0.25 / 0.125 / 0.0625`。

诊断显示多个固定组分内部有 5 个 distinct `rho`，尤其：

- `3c_B_Ch:16,4c_B_Ch:8,4c_Cd_1_Ch:4`: `n=49`，`distinct_rho=5`，`has_High=True`。
- `3c_B_Ch:16,4c_B_Ch:4,4c_Cd_1_Ch:8`: `n=49`，`distinct_rho=5`，`has_High=True`。

因此 `fjh` 可以真正检验 composition-controlled symmetry-energy trend。

## Current Selection Strategy

1. 先按 topology 枚举所有 symmetry-inequivalent colorings。
2. 按 composition，即 ligand count fingerprint，分组。
3. 只保留组分内 distinct `rho` 足够多的组。
4. `ffc`: 成本低，穷举所有组分内 `rho >= 2` 的组分。
5. `fjh`: 主验证证据，成本高，挑 `rho` 最丰富的组分，按每个 distinct `rho` 分层抽样。当前建议用 2 个 composition、每个 `rho` 最多取 4 个，以提高最低能点和 Spearman 趋势的可信度。
6. 当前论文论点不是“`rho` 越高能量越低”的单调关系，而是“低能候选多半保留在中/高对称区域”，即 symmetry-prioritized enumeration 可作为 pre-screening heuristic。`ffc` 可 pooled 作描述性图；composition-controlled 分析用于避免组分混杂。

这个设计比直接把所有结构混在一起更科学，因为能量强烈受局部组分和配体摆放影响。

## CP2K / DFT Method

服务器：

```text
huangbo@Lab207
112 cores / 128 GB RAM
No GPU
Offline
Shared by multiple users
```

CP2K binary:

```text
/home/huangbo/cp2k/cp2k-2025.2-Linux-gnu-x86_64.ssmp
```

CP2K data directory:

```text
/home/huangbo/cp2k/cp2k-2025.2/data
```

当前使用：

- Pure OpenMP, no MPI。
- Batch-level concurrency: multiple independent CP2K single-point jobs can run in parallel by topology.
- `BASIS_MOLOPT_UZH`。
- `GTH_POTENTIALS`。
- `dftd3.dat`。
- PBE+D3。
- `DZVP-MOLOPT-PBE-GTH` basis。
- GTH pseudopotentials。
- cutoff / relative cutoff: `500 / 60 Ry`。

由于 transition metals + conjugated linkers 小能隙，OT SCF 不收敛，当前已改为：

- Diagonalization。
- Fermi-Dirac smearing at `300 K`。
- `ADDED_MOS 30`。
- Broyden mixing。
- `EPS_SCF = 5e-5`。
- `MAX_SCF = 300`。
- `IGNORE_CONVERGENCE_FAILURE = True`。

并发配置位于 `dft_validation/config.py`：

```python
CP2K_OMP = 16
DFT_PARALLEL = True
DFT_MAX_WORKERS = {"ffc": 4, "fjh": 2}
DFT_JOB_MEM_GB = {"ffc": 24, "fjh": 55}
```

在 112 cores / 128 GB RAM 的服务器上，这意味着 `ffc` 最多约 64 个活跃 OpenMP 线程，`fjh` 最多约 32 个活跃 OpenMP 线程。`fjh` 保守并发是因为 760 原子对角化非常吃内存。

当前只做 single-point `ENERGY`，不做 relaxation。论文口径是：DFT 用于评估 enumeration-time pre-screening score，而非最终 thermodynamic stability；固定刚性组装初始几何可避免 relaxation path 混杂。报告能量应称为 geometry-consistent screening scores, not relaxed formation energies。

## Server Runbook

启动或继续 DFT：

```bash
cd ~/mof
export PYTHONPATH=$PWD:$PYTHONPATH
nohup python -u dft_validation/run_dft_validation.py --resume > run_all.log 2>&1 &
tail -f run_all.log
```

查看当前 CP2K 输出路径：

```bash
ps aux | grep '[c]p2k'
```

然后 tail 其中 `-o` 后面的 `cp2k_run.out`。

注意事项：

- `Ctrl+C` 停止 `tail` 不会停止 nohup Python。
- 要停任务时，先 kill 外层 Python，再确认 CP2K 没残留。
- 不要在 CP2K 运行中删除 `runs/` 或 `relaxed/`。
- 如果要部署新的并发代码，先停止旧的串行 nohup Python 和它启动的 CP2K，避免同一结构重复跑。

## Paper Direction

`paper/main.tex` 里的 Symmetry-Stability Correlation 和 DFT Energy Validation 已改成更稳的表述：

- 不再声称 symmetry 一定导致更低能。
- 改为 composition-controlled qualitative validation。
- 报告低能候选是否落在按固定 `rho` 区间定义的中/高对称区域；composition-controlled 子集可额外报告 Spearman correlation 和 maximum-`rho` 检查。
- 论点是 symmetry-prioritized enumeration 是 a priori pre-screening heuristic，而不是 stability guarantee。
- 当前论文图件使用并列 pooled 图：`paper/ffc_pooled_symmetry_energy.png` 和 `paper/fjh_pooled_symmetry_energy.png`。原 `fig:symmetry_stability` composition-controlled 双图已从论文主文删除，避免重复叙事。
- pooled 图是主文中的描述性证据：实心点展示所有已完成候选，空心点展示 fixed composition 内最低能候选。论文结论应基于这些 composition-wise minima 是否落在 medium/high symmetry 区间，而不是基于全局 pooled slope 或“rho-energy 单调下降”。
- 已跑完的 CP2K 不应浪费：使用 `python -u dft_validation/run_dft_validation.py --resume --harvest-only --summary` 可从已有 `runs/*/cp2k_run.out` 解析能量，更新 JSON 并重新分析，不会启动新 CP2K。
- pooled 图均按统一阈值绘制：`0 < rho <= 0.25` Low, `0.25 < rho <= 0.75` Medium, `0.75 < rho <= 1` High。空心圆标记每个 fixed composition 内最低能候选。`fjh` 当前图只有已完成计算子集，应在服务器完成更多 fjh 点后更新。

继续工作时必须保持论文与代码一致：

- PBE+D3。
- `500 / 60 Ry`。
- Single-point energy。
- Composition-controlled stratified sampling。
