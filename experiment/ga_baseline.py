import os
import sys
import json
import random
import time
import csv
from deap import base, creator, tools

# [配置] 解决 macOS 下 matplotlib 报错并设置为后台绘图
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# [新增导入] CHGNet 和 pymatgen
from chgnet.model import CHGNet
from pymatgen.core import Structure

# 路径修正：确保能导入项目中的自定义模块
current_dir = os.path.dirname(os.path.abspath(__file__))
if os.path.exists(os.path.join(current_dir, "../data")):
    project_root = current_dir
else:
    project_root = os.path.dirname(current_dir)
sys.path.append(project_root)

try:
    from symmetry import SymmetryEngine
    from enumerator import MOFEnumerator
    from assembler import StructureAssembler  # [新增导入] 引入组装器
except ImportError:
    print("Warning: Custom modules 'symmetry', 'enumerator' or 'assembler' not found.")


def load_graph(path):
    with open(path, 'r') as f: return json.load(f)


# ================= 辅助函数：事后正则化 =================
def _apply_perm(coloring, perm):
    n = len(coloring)
    new_c = [""] * n
    for i in range(n):
        new_c[perm[i]] = coloring[i]
    return tuple(new_c)


def get_canonical(coloring, perms):
    """计算构型的正则代表元，用于严格消除同构冗余"""
    canonical = coloring
    for p in perms:
        candidate = _apply_perm(coloring, p)
        if candidate < canonical:
            canonical = candidate
    return canonical


# ================= 核心 GA 评估函数 (替换为 CHGNet) =================
def evaluate_chgnet_energy(individual, graph_data, atom_choices, assembler, chgnet_model):
    """使用 CHGNet 计算重构晶体的单原子能量"""
    coloring = [atom_choices[idx] for idx in individual]
    structure_map = {str(i): color for i, color in enumerate(coloring)}

    # 生成一个随机的临时文件名，防止覆盖冲突
    temp_cif_path = os.path.join(assembler.output_dir, f"temp_ga_{random.randint(100000, 999999)}.cif")

    try:
        # 1. 组装 3D 结构
        assembler.assemble_and_save(
            filename=temp_cif_path,
            structure_map=structure_map,
            graph_data=graph_data
        )

        # 2. 加载结构并计算能量
        struct = Structure.from_file(temp_cif_path)
        prediction = chgnet_model.predict_structure(struct)
        energy_per_atom = prediction['e']

        # 3. 评估成功，清理临时文件
        if os.path.exists(temp_cif_path): os.remove(temp_cif_path)
        log_file = temp_cif_path.replace(".cif", ".log")
        if os.path.exists(log_file): os.remove(log_file)

        return (energy_per_atom,)

    except Exception as e:
        # 拓扑上存在严重的几何重叠或物理不合理导致组装/预测失败
        if os.path.exists(temp_cif_path): os.remove(temp_cif_path)
        # 给予一个极高的能量惩罚，让 GA 在选择阶段自动淘汰该个体
        return (9999.0,)


def main():
    # ================= 0. 初始化深度学习模型与组装器 =================
    print(">>> Loading CHGNet Universal Potential...")
    chgnet_model = CHGNet.load()

    ga_temp_dir = os.path.join(project_root, "../results", "ga_temp_cifs")
    os.makedirs(ga_temp_dir, exist_ok=True)
    assembler = StructureAssembler(output_dir=ga_temp_dir)

    # ================= 1. 实验配置区 =================
    TOPOLOGY = "sod"
    # ATOM_CHOICES = ["6c_CdZn_1_Ch", "8c_Co_1_Ch", "0_mn1_new"]
    ATOM_CHOICES = ["8c_Co_1_Ch", "0_mn1_new"]

    NUM_GENERATIONS = 150
    POPULATION_SIZE = 300
    ELITISM_SIZE = 10  # 精英保留数量
    CXPB = 0.7  # 交叉概率
    MUTPB = 0.2  # 变异概率

    atom_count = len(ATOM_CHOICES)
    # 估算的总评估次数（仅用于输出文件命名匹配）
    TOTAL_EVALS_ESTIMATE = POPULATION_SIZE * NUM_GENERATIONS

    graph_path = os.path.join(project_root, "../data", TOPOLOGY, f"{TOPOLOGY}_graph.json")
    if not os.path.exists(graph_path):
        print(f"Error: Graph file not found at {graph_path}")
        return

    print(f"--- Loading Topology: {TOPOLOGY} ---")
    graph_data = load_graph(graph_path)
    num_nodes = graph_data["n_nodes"]

    # ================= 2. 获取对称性与 Ground Truth =================
    print("\n--- Step 1: Pre-computing Symmetry and Ground Truth ---")
    sym_engine = SymmetryEngine(graph_data, verbose=False)
    perms = sym_engine.compute_automorphism_group()
    print(f"  Automorphism Group Size |G|: {len(perms)}")

    print("  Computing Ground Truth (Burnside/Stratified Enumerator)...")
    enum = MOFEnumerator(perms, num_nodes, verbose=False)
    all_unique_structures = enum.enumerate(ATOM_CHOICES)

    # 获取目标高对称性结构 (Stabilizer > 1)
    high_sym_canonicals = set()
    for res in all_unique_structures:
        canon = tuple(res['coloring'])
        if res['stabilizer_size'] > 1:
            high_sym_canonicals.add(canon)

    total_high_sym_targets = len(high_sym_canonicals)
    print(f"   -> Ground Truth: {len(all_unique_structures)} unique orbits.")
    print(f"   -> High-Sym Targets (Stab > 1): {total_high_sym_targets}")

    if total_high_sym_targets == 0:
        print("   -> Warning: No High-Symmetry targets found in Ground Truth.")

    # ================= 3. DEAP GA 环境设置 =================
    print("\n--- Step 2: Initializing Advanced Genetic Algorithm ---")
    if not hasattr(creator, "FitnessMin"):
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    if not hasattr(creator, "Individual"):
        creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()
    toolbox.register("attr_atom", random.randint, 0, len(ATOM_CHOICES) - 1)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_atom, n=num_nodes)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # 注册修改后的 CHGNet 评估函数
    toolbox.register("evaluate", evaluate_chgnet_energy,
                     graph_data=graph_data,
                     atom_choices=ATOM_CHOICES,
                     assembler=assembler,
                     chgnet_model=chgnet_model)

    toolbox.register("mate", tools.cxUniform, indpb=0.5)
    toolbox.register("mutate", tools.mutUniformInt, low=0, up=len(ATOM_CHOICES) - 1, indpb=0.1)
    toolbox.register("select", tools.selTournament, tournsize=3)

    # ================= 4. 带严密监控的进化循环 =================
    print(f"  Running GA: Pop={POPULATION_SIZE}, Gens={NUM_GENERATIONS}, Elitism={ELITISM_SIZE}...")
    pop = toolbox.population(n=POPULATION_SIZE)

    seen_canonicals = set()
    unique_ground_states = set()

    # [新增] 能量缓存字典：避免对同构重复进行费时的 3D 组装和预测
    energy_cache = {}

    total_evals = 0
    redundant_evals = 0
    global_min_energy = float('inf')

    csv_data = []

    # 记录作图所需数组
    history_evals = []
    history_high_cov = []

    start_time = time.time()

    def process_evaluation(ind):
        """核心评估监控逻辑，集成事后正则化与同构缓存机制"""
        nonlocal total_evals, redundant_evals, global_min_energy, unique_ground_states
        total_evals += 1

        coloring = tuple(ATOM_CHOICES[idx] for idx in ind)
        canon = get_canonical(coloring, perms)

        # 【核心逻辑】：拦截同构冗余，直接命中缓存，大幅节约算力
        if canon in energy_cache:
            redundant_evals += 1
            ind.fitness.values = (energy_cache[canon],)  # 命中缓存
        else:
            # 只有遇到物理上全新的构型，才执行深度学习预测
            ind.fitness.values = toolbox.evaluate(ind)
            energy = ind.fitness.values[0]

            # 记录到缓存
            energy_cache[canon] = energy
            seen_canonicals.add(canon)

            if energy < global_min_energy:
                global_min_energy = energy
                unique_ground_states = {canon}
            elif energy == global_min_energy:
                unique_ground_states.add(canon)

    # 4.1 评估初始种群
    for ind in pop:
        process_evaluation(ind)

    # 4.2 初始种群指标记录
    current_high_sym_count = sum(1 for c in seen_canonicals if c in high_sym_canonicals)
    current_coverage_pct = (current_high_sym_count / total_high_sym_targets * 100) if total_high_sym_targets > 0 else 0

    csv_data.append({
        "Generation": 0,
        "Total_Evaluations": total_evals,
        "Unique_Explored": len(seen_canonicals),
        "Redundant_Evaluations": redundant_evals,
        "Redundancy_Rate": redundant_evals / total_evals if total_evals > 0 else 0,
        "Best_Energy": global_min_energy,
        "High_Sym_Coverage_Pct": current_coverage_pct
    })
    history_evals.append(total_evals)
    history_high_cov.append(current_coverage_pct)

    # 4.3 正式迭代
    for gen in range(1, NUM_GENERATIONS + 1):
        # 精英保留
        elites = tools.selBest(pop, ELITISM_SIZE)
        offspring = list(map(toolbox.clone, elites))

        # 锦标赛选择
        mating_pool = toolbox.select(pop, POPULATION_SIZE - ELITISM_SIZE)
        mating_pool = list(map(toolbox.clone, mating_pool))

        # 交叉
        for child1, child2 in zip(mating_pool[::2], mating_pool[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        # 变异
        for mutant in mating_pool:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        offspring.extend(mating_pool)

        # 评估无效个体（新生成的个体）
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        for ind in invalid_ind:
            process_evaluation(ind)

        pop[:] = offspring

        # 当前代指标统计
        current_high_sym_count = sum(1 for c in seen_canonicals if c in high_sym_canonicals)
        current_coverage_pct = (
                    current_high_sym_count / total_high_sym_targets * 100) if total_high_sym_targets > 0 else 0
        redundancy_rate = redundant_evals / total_evals if total_evals > 0 else 0

        csv_data.append({
            "Generation": gen,
            "Total_Evaluations": total_evals,
            "Unique_Explored": len(seen_canonicals),
            "Redundant_Evaluations": redundant_evals,
            "Redundancy_Rate": redundancy_rate,
            "Best_Energy": global_min_energy,
            "High_Sym_Coverage_Pct": current_coverage_pct
        })
        history_evals.append(total_evals)
        history_high_cov.append(current_coverage_pct)

        if gen % 10 == 0 or gen == NUM_GENERATIONS:
            print(f"  Gen {gen}: Evals={total_evals} | Best Energy={global_min_energy:.4f} | "
                  f"High-Sym Cov={current_coverage_pct:.2f}% | Redundancy={redundancy_rate:.1%}")

    ga_time = time.time() - start_time
    print(f"\n  GA finished in {ga_time:.2f} seconds.")

    # ================= 5. 导出 CSV 与 绘图 =================
    output_dir = os.path.join(project_root, "../results", "ga_baseline")
    os.makedirs(output_dir, exist_ok=True)

    # 命名风格: experiment_data_{TOPOLOGY}_{TOTAL_EVALS}_{atom_count}
    filename_base = f"experiment_data_{TOPOLOGY}_{TOTAL_EVALS_ESTIMATE}_{atom_count}"
    csv_path = os.path.join(output_dir, f"{filename_base}.csv")

    fieldnames = ["Generation", "Total_Evaluations", "Unique_Explored",
                  "Redundant_Evaluations", "Redundancy_Rate",
                  "Best_Energy", "High_Sym_Coverage_Pct"]

    print(f"\n>>> Saving data to {csv_path}...")
    with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in csv_data:
            # 格式化输出保证整洁
            row_copy = row.copy()
            row_copy['Redundancy_Rate'] = f"{row['Redundancy_Rate']:.4f}"
            row_copy['High_Sym_Coverage_Pct'] = f"{row['High_Sym_Coverage_Pct']:.4f}"
            writer.writerow(row_copy)

    print(">>> Plotting...")
    plt.figure(figsize=(4, 4))
    plt.subplot(1, 1, 1)

    # 绘制高对称覆盖率曲线
    plt.plot(history_evals, history_high_cov, color='g', label='Genetic Algorithm')

    plt.xlabel('Number of Evaluations')
    plt.ylabel('Symmetric Target Coverage (%)')
    plt.title(f'{TOPOLOGY}, {atom_count} atoms')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    img_path = os.path.join(output_dir, f"{filename_base}.png")
    plt.savefig(img_path, dpi=300)
    print(f"   -> Plot Saved to {img_path}")
    print(f"   -> Final Redundancy Rate: {csv_data[-1]['Redundancy_Rate']:.1%}")


if __name__ == "__main__":
    main()