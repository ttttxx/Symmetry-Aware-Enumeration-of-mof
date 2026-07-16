import os
import sys
import json
import random
import time
import csv
import math
from deap import base, creator, tools

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from chgnet.model import CHGNet
from pymatgen.core import Structure

current_dir = os.path.dirname(os.path.abspath(__file__))
if os.path.exists(os.path.join(current_dir, "../data")):
    project_root = current_dir
else:
    project_root = os.path.dirname(current_dir)
sys.path.append(project_root)

try:
    from symmetry import SymmetryEngine
    from enumerator import MOFEnumerator
    from assembler import StructureAssembler
    from main import calculate_topology_requirements
except ImportError:
    print("Warning: Custom modules 'symmetry', 'enumerator' or 'assembler' not found.")

def load_graph(path):
    with open(path, 'r') as f: return json.load(f)

def _apply_perm(coloring, perm):
    n = len(coloring)
    new_c = [""] * n
    for i in range(n):
        new_c[perm[i]] = coloring[i]
    return tuple(new_c)

def get_canonical(coloring, perms):
    canonical = coloring
    for p in perms:
        candidate = _apply_perm(coloring, p)
        if candidate < canonical:
            canonical = candidate
    return canonical

def evaluate_chgnet_energy(individual, graph_data, atom_choices, assembler, chgnet_model):
    coloring = [atom_choices[idx] for idx in individual]
    structure_map = {str(i): color for i, color in enumerate(coloring)}

    temp_cif_path = os.path.join(assembler.output_dir, f"temp_ga_{random.randint(100000, 999999)}.cif")

    try:

        assembler.assemble_and_save(
            filename=temp_cif_path,
            structure_map=structure_map,
            graph_data=graph_data
        )

        struct = Structure.from_file(temp_cif_path)
        prediction = chgnet_model.predict_structure(struct)
        energy_per_atom = prediction['e']

        if os.path.exists(temp_cif_path): os.remove(temp_cif_path)
        log_file = temp_cif_path.replace(".cif", ".log")
        if os.path.exists(log_file): os.remove(log_file)

        return (energy_per_atom,)

    except Exception as e:

        if os.path.exists(temp_cif_path): os.remove(temp_cif_path)

        return (9999.0,)

def compute_ground_truth_min_energy(all_unique_structures, graph_data, assembler,
                                    chgnet_model, max_eval=2000):
    n = len(all_unique_structures)
    if n > max_eval:
        print(f"  [GT-Energy] Skipped: {n} structures exceed max_eval={max_eval}; "
              f"reference line will be unavailable for this case.")
        return (float('nan'), None, n)

    best_e = float('inf')
    best_stab = None
    for k, res in enumerate(all_unique_structures):
        coloring = res['coloring']
        structure_map = {str(i): color for i, color in enumerate(coloring)}
        temp_cif = os.path.join(assembler.output_dir, f"temp_gt_{k}.cif")
        try:
            assembler.assemble_and_save(filename=temp_cif,
                                        structure_map=structure_map,
                                        graph_data=graph_data)
            struct = Structure.from_file(temp_cif)
            e = chgnet_model.predict_structure(struct)['e']
        except Exception:
            e = 9999.0
        finally:
            if os.path.exists(temp_cif):
                os.remove(temp_cif)
            logf = temp_cif.replace(".cif", ".log")
            if os.path.exists(logf):
                os.remove(logf)

        if e < best_e:
            best_e = e
            best_stab = res.get('stabilizer_size', None)

    return (best_e, best_stab, n)

def main():

    print(">>> Loading CHGNet Universal Potential...")
    chgnet_model = CHGNet.load()

    ga_temp_dir = os.path.join(project_root, "../results", "ga_temp_cifs")
    os.makedirs(ga_temp_dir, exist_ok=True)
    assembler = StructureAssembler(output_dir=ga_temp_dir)

    test_cases = [

        # # 3c×1, 4c×2 → naive=1^16*2^12=4096
        # # 3c×2, 4c×1 → naive=2^16*1^12=65536
        # # 3c×1, 4c×3 → naive=1^16*3^12=531441
        ("fjh", {3: ["3c_B_Ch", "3c_triazine_Ch"], 4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"]}),

    ]
    TOPOLOGY, USER_SELECTED_LIGANDS = test_cases[0]
    ATOM_CHOICES = sorted(set(lig for ligs in USER_SELECTED_LIGANDS.values() for lig in ligs))

    NUM_GENERATIONS = 150
    POPULATION_SIZE = 300
    ELITISM_SIZE = 10
    CXPB = 0.7
    MUTPB = 0.2

    atom_count = len(ATOM_CHOICES)

    TOTAL_EVALS_ESTIMATE = POPULATION_SIZE * NUM_GENERATIONS

    graph_path = os.path.join(project_root, "../data", TOPOLOGY, f"{TOPOLOGY}_graph.json")
    if not os.path.exists(graph_path):
        print(f"Error: Graph file not found at {graph_path}")
        return

    print(f"--- Loading Topology: {TOPOLOGY} ---")
    graph_data = load_graph(graph_path)
    num_nodes = graph_data["n_nodes"]
    node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]
    topology_requirements = calculate_topology_requirements(graph_data)

    print("\n--- Step 1: Pre-computing Symmetry and Ground Truth ---")
    sym_engine = SymmetryEngine(graph_data, verbose=False)
    perms = sym_engine.compute_automorphism_group()
    print(f"  Automorphism Group Size |G|: {len(perms)}")

    print("  Computing Ground Truth (Burnside/Stratified Enumerator)...")
    enum = MOFEnumerator(perms, num_nodes, verbose=False,
                         topology_requirements=topology_requirements,
                         user_selected_ligands=USER_SELECTED_LIGANDS,
                         node_degrees=node_degrees)
    all_unique_structures = enum.enumerate(ATOM_CHOICES)

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

    print("  Computing exhaustive global-minimum energy (reference for energy-convergence panel)...")
    gt_min_energy, gt_min_stab, gt_n = compute_ground_truth_min_energy(
        all_unique_structures, graph_data, assembler, chgnet_model
    )
    if not math.isnan(gt_min_energy):
        print(f"   -> Enumerated global minimum = {gt_min_energy:.4f} eV/atom "
              f"(|Stab|={gt_min_stab}, over {gt_n} structures)")

    print("\n--- Step 2: Initializing Advanced Genetic Algorithm ---")
    if not hasattr(creator, "FitnessMin"):
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    if not hasattr(creator, "Individual"):
        creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()
    toolbox.register("attr_atom", random.randint, 0, len(ATOM_CHOICES) - 1)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_atom, n=num_nodes)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("evaluate", evaluate_chgnet_energy,
                     graph_data=graph_data,
                     atom_choices=ATOM_CHOICES,
                     assembler=assembler,
                     chgnet_model=chgnet_model)

    toolbox.register("mate", tools.cxUniform, indpb=0.5)
    toolbox.register("mutate", tools.mutUniformInt, low=0, up=len(ATOM_CHOICES) - 1, indpb=0.1)
    toolbox.register("select", tools.selTournament, tournsize=3)

    print(f"  Running GA: Pop={POPULATION_SIZE}, Gens={NUM_GENERATIONS}, Elitism={ELITISM_SIZE}...")
    pop = toolbox.population(n=POPULATION_SIZE)

    seen_canonicals = set()
    unique_ground_states = set()

    energy_cache = {}

    total_evals = 0
    redundant_evals = 0
    global_min_energy = float('inf')

    csv_data = []

    history_evals = []
    history_high_cov = []

    start_time = time.time()

    def process_evaluation(ind):
        nonlocal total_evals, redundant_evals, global_min_energy, unique_ground_states
        total_evals += 1

        coloring = tuple(ATOM_CHOICES[idx] for idx in ind)
        canon = get_canonical(coloring, perms)

        if canon in energy_cache:
            redundant_evals += 1
            ind.fitness.values = (energy_cache[canon],)
        else:

            ind.fitness.values = toolbox.evaluate(ind)
            energy = ind.fitness.values[0]

            energy_cache[canon] = energy
            seen_canonicals.add(canon)

            if energy < global_min_energy:
                global_min_energy = energy
                unique_ground_states = {canon}
            elif energy == global_min_energy:
                unique_ground_states.add(canon)

    for ind in pop:
        process_evaluation(ind)

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

    for gen in range(1, NUM_GENERATIONS + 1):

        elites = tools.selBest(pop, ELITISM_SIZE)
        offspring = list(map(toolbox.clone, elites))

        mating_pool = toolbox.select(pop, POPULATION_SIZE - ELITISM_SIZE)
        mating_pool = list(map(toolbox.clone, mating_pool))

        for child1, child2 in zip(mating_pool[::2], mating_pool[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        for mutant in mating_pool:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        offspring.extend(mating_pool)

        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        for ind in invalid_ind:
            process_evaluation(ind)

        pop[:] = offspring

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

    output_dir = os.path.join(project_root, "../results", "ga_baseline")
    os.makedirs(output_dir, exist_ok=True)

    filename_base = f"experiment_data_{TOPOLOGY}_{TOTAL_EVALS_ESTIMATE}_{atom_count}"
    csv_path = os.path.join(output_dir, f"{filename_base}.csv")

    fieldnames = ["Generation", "Total_Evaluations", "Unique_Explored",
                  "Redundant_Evaluations", "Redundancy_Rate",
                  "Best_Energy", "High_Sym_Coverage_Pct",
                  "GT_Min_Energy", "GT_Min_Stab"]

    gt_energy_str = "" if math.isnan(gt_min_energy) else f"{gt_min_energy:.6f}"
    gt_stab_str = "" if gt_min_stab is None else str(gt_min_stab)

    print(f"\n>>> Saving data to {csv_path}...")
    with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in csv_data:

            row_copy = row.copy()
            row_copy['Redundancy_Rate'] = f"{row['Redundancy_Rate']:.4f}"
            row_copy['High_Sym_Coverage_Pct'] = f"{row['High_Sym_Coverage_Pct']:.4f}"
            row_copy['GT_Min_Energy'] = gt_energy_str
            row_copy['GT_Min_Stab'] = gt_stab_str
            writer.writerow(row_copy)

    print(">>> Plotting...")
    plt.figure(figsize=(4, 4))
    plt.subplot(1, 1, 1)

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
    if not math.isnan(gt_min_energy):
        print(f"   -> GA best = {global_min_energy:.4f} | Enumerated global min = {gt_min_energy:.4f} "
              f"| gap = {global_min_energy - gt_min_energy:.4f} eV/atom")

if __name__ == "__main__":
    main()