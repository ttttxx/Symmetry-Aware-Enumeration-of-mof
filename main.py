# mof/main.py
import os
import json
import random
import numpy as np
from datetime import datetime
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator
from verification import BridgeVerifier
from assembler import StructureAssembler


def enrich_graph_data_with_topology(graph_data):
    """
    [终极鲁棒版]：直接对接 parse_cif 的完美拓扑，废弃有缺陷的坐标系去重猜测。
    利用真实的晶格平移 sym 计算 Cycle Basis，通杀任何扭曲或多重边的 MOF。
    """
    nodes = graph_data['nodes']
    unique_edges = graph_data['edges']
    V = len(nodes)
    E = len(unique_edges)
    node_idx = {str(n.get('id', i)): i for i, n in enumerate(nodes)}

    print(f"  [Topology] Using {E} robust edges from parse_cif for {V} nodes.")

    omega = np.zeros((E, 3))
    for k, edge in enumerate(unique_edges):
        u, v = node_idx[str(edge['source'])], node_idx[str(edge['target'])]
        f_u, f_v = np.array(nodes[u]['frac']), np.array(nodes[v]['frac'])

        sym = edge.get('sym', '.')
        if sym == '.':
            shift = np.array([0, 0, 0])
        else:
            s_str = sym.split('_')[1]
            shift = np.array([int(s_str[0]) - 5, int(s_str[1]) - 5, int(s_str[2]) - 5])

        omega[k] = f_v + shift - f_u

    parent = list(range(V))

    def find(i):
        if parent[i] == i: return i
        parent[i] = find(parent[i])
        return parent[i]

    tree_edges, chords = [], []
    for k, edge in enumerate(unique_edges):
        u, v = node_idx[str(edge['source'])], node_idx[str(edge['target'])]
        root_u, root_v = find(u), find(v)
        if root_u != root_v:
            parent[root_u] = root_v
            tree_edges.append(k)
        else:
            chords.append(k)

    C = np.zeros((len(chords), E))
    T_adj = {i: [] for i in range(V)}
    for k in tree_edges:
        u, v = node_idx[str(unique_edges[k]['source'])], node_idx[str(unique_edges[k]['target'])]
        T_adj[u].append((v, k, 1))
        T_adj[v].append((u, k, -1))

    def find_path(start, end):
        queue = [(start, [])]
        visited = {start}
        while queue:
            curr, path = queue.pop(0)
            if curr == end: return path
            for nxt, idx, d in T_adj[curr]:
                if nxt not in visited:
                    visited.add(nxt)
                    queue.append((nxt, path + [(idx, d)]))
        return None

    alpha_matrix = np.zeros((len(chords), 3))
    for i, c_idx in enumerate(chords):
        u, v = node_idx[str(unique_edges[c_idx]['source'])], node_idx[str(unique_edges[c_idx]['target'])]
        C[i, c_idx] = 1
        path = find_path(v, u)
        if path:
            for e_idx, d in path: C[i, e_idx] = d
        alpha_matrix[i] = np.dot(C[i], omega)

    M_all = np.zeros((V - 1, E))
    for k, edge in enumerate(unique_edges):
        u, v = node_idx[str(edge['source'])], node_idx[str(edge['target'])]
        if u < V - 1: M_all[u, k] = -1
        if v < V - 1: M_all[v, k] = 1

    graph_data['Bstar'] = np.vstack([C, M_all]).tolist()
    graph_data['alpha_matrix'] = np.round(alpha_matrix).astype(int).tolist()
    graph_data['n_edges'] = E
    return graph_data


def load_graph(path):
    with open(path, 'r') as f: return json.load(f)


def calculate_topology_requirements(graph_data):
    """
    Analyze degree distribution in the topology graph.
    Returns dict like {3: 4, 4: 3} — need 4 three-c nodes and 3 four-c nodes.
    """
    degree_counts = {}
    for node in graph_data['nodes']:
        d = node.get('degree', 0)
        degree_counts[d] = degree_counts.get(d, 0) + 1
    return degree_counts


def generate_output_filename(prefix, atom_names):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    atoms_str = "".join(sorted(atom_names))
    filename = f"{prefix}_{atoms_str}_{timestamp}.json"
    return os.path.join("results", filename)


# ══════════════════════════════════════════════════════════════════════════════
# 对称性分层抽样
# ══════════════════════════════════════════════════════════════════════════════

def sample_by_symmetry(results, group_size,
                        n_per_cat=20,
                        high_rho=0.75, low_rho=0.25,
                        seed=42):
    """
    将枚举结果按归一化对称指数 ρ = |Stab(c)| / |Γ| 分为三类，
    每类随机抽取最多 n_per_cat 个结构。

    分类规则（与 run_analysis_plot_CHGnet.py 保持一致）：
      High : ρ ∈ [high_rho, 1.0]   — 高度有序，对称性最高
      Med  : ρ ∈ (low_rho, high_rho) — 中等对称
      Low  : ρ ∈ [0,  low_rho]     — 低对称 / 非晶态

    参数
    ----
    results    : enumerator 返回的结构列表（含 'stabilizer_size', 'coloring'）
    group_size : 自同构群阶数 |Γ|
    n_per_cat  : 每类最多抽取的结构数量
    high_rho   : High/Med 分界阈值
    low_rho    : Med/Low  分界阈值
    seed       : 随机种子，保证可重现

    返回
    ----
    list of (res_data, category_str, rho_value)
    """
    rng = random.Random(seed)

    pools = {'High': [], 'Med': [], 'Low': []}
    for r in results:
        rho = r['stabilizer_size'] / group_size
        if rho >= high_rho:
            pools['High'].append(r)
        elif rho > low_rho:
            pools['Med'].append(r)
        else:
            pools['Low'].append(r)

    sampled = []
    for cat in ['High', 'Med', 'Low']:
        pool = pools[cat]
        chosen = pool if len(pool) <= n_per_cat else rng.sample(pool, n_per_cat)
        print(f"    [{cat:4s}] pool={len(pool):5d}  sampled={len(chosen)}")
        for r in chosen:
            rho = r['stabilizer_size'] / group_size
            sampled.append((r, cat, rho))

    total = len(sampled)
    print(f"    → 共抽取 {total} 个结构"
          f"（High={len([s for s in sampled if s[1]=='High'])}"
          f" Med={len([s for s in sampled if s[1]=='Med'])}"
          f" Low={len([s for s in sampled if s[1]=='Low'])}）")
    return sampled


def run_enumeration(data_prefix, user_selected_ligands, forced_edge,
                    max_cifs=5, max_results=10_000,
                    sample_mode=False, n_per_cat=20,
                    high_rho=0.75, low_rho=0.25):
    """
    Parameterized entry point for the Streamlit frontend.
    user_selected_ligands: dict[int, list[str]] — e.g. {3: ["3c_B_Ch"], 4: ["4c_Cu_1_Ch"]}

    sample_mode : False（默认）→ 按 stabilizer_size 排序，输出前 max_cifs 个（原行为）
                  True         → 按 ρ 阈值分层，每类输出最多 n_per_cat 个 CIF
    n_per_cat   : sample_mode=True 时每类的抽样上限
    high_rho    : High/Med 分界（默认 0.75）
    low_rho     : Med/Low  分界（默认 0.25）

    Returns: (cif_paths: list[str], json_path: str, output_data: dict)
    """
    graph_path = os.path.join("data", data_prefix, f"{data_prefix}_graph.json")
    if not os.path.exists(graph_path):
        raise FileNotFoundError(f"Graph file not found at {graph_path}")

    print(f"--- Loading Topology: {data_prefix} ---")
    graph_data = load_graph(graph_path)
    graph_data = enrich_graph_data_with_topology(graph_data)
    num_nodes = graph_data["n_nodes"]
    node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]

    # ---- 自动推导拓扑需求 ----
    topology_requirements = calculate_topology_requirements(graph_data)
    print(f"  [Topology] Requirements: {topology_requirements}")

    # 展平得到全部配体名 (用于输出文件名与 enumerator 的 flat list)
    all_ligands = sorted(set(
        lig for ligs in user_selected_ligands.values() for lig in ligs
    ))

    print("\n--- Step 1: Symmetry Analysis ---")
    sym_engine = SymmetryEngine(graph_data, verbose=True)
    perms = sym_engine.compute_automorphism_group()
    group_size = len(perms)
    print(f"  Automorphism Group Size |G|: {group_size}")

    print("\n--- Step 2: Structure Enumeration ---")
    enumerator = MOFEnumerator(perms, num_nodes, verbose=True,
                               topology_requirements=topology_requirements,
                               user_selected_ligands=user_selected_ligands,
                               node_degrees=node_degrees)
    results = enumerator.enumerate(all_ligands)
    print(f"  Found {len(results)} unique structures (early-pruned by degree & count constraints).")
    if len(results) > max_results:
        print(f"  Limiting to {max_results} results.")
        results = results[:max_results]

    print("\n--- Step 3: Verification ---")
    verifier = BridgeVerifier(group_size)
    is_verified, stats = verifier.verify(results, num_nodes, all_ligands, None,
                                          topology_requirements=topology_requirements,
                                          user_selected_ligands=user_selected_ligands)

    print("\n--- Step 4: Assembler (Generate CIFs) ---")
    atoms_str = "".join(sorted(all_ligands))
    cif_subdir = f"{data_prefix}_{atoms_str}_cifs"
    cif_output_path = os.path.join("results", cif_subdir)
    os.makedirs(cif_output_path, exist_ok=True)

    assembler = StructureAssembler(output_dir=cif_output_path)
    cif_paths = []

    if sample_mode:
        # ── 分层抽样模式 ────────────────────────────────────────────────────
        print(f"  [采样模式] ρ 阈值: High≥{high_rho} / Med ({low_rho},{high_rho}) / Low≤{low_rho}")
        print(f"  [采样模式] 每类最多 {n_per_cat} 个 CIF …")
        sampled = sample_by_symmetry(results, group_size,
                                      n_per_cat=n_per_cat,
                                      high_rho=high_rho, low_rho=low_rho)
        for i, (res_data, cat, rho) in enumerate(sampled):
            stab = res_data['stabilizer_size']
            fname = os.path.join(
                cif_output_path,
                f"{data_prefix}_{atoms_str}_{cat}_sym{stab}_rho{rho:.3f}_{i:03d}.cif")
            smap = {str(j): elem for j, elem in enumerate(res_data['coloring'])}
            try:
                assembler.assemble_and_save(fname, smap, graph_data, forced_edge=forced_edge)
                cif_paths.append(fname)
            except Exception as e:
                print(f"  [Error] {os.path.basename(fname)}: {e}")
    else:
        # ── 原始 Top-K 模式 ─────────────────────────────────────────────────
        results.sort(key=lambda x: x['stabilizer_size'], reverse=True)
        print(f"  Generating CIFs for Top-{min(len(results), max_cifs)} structures...")
        for i, res_data in enumerate(results[:max_cifs]):
            stab_size = res_data['stabilizer_size']
            fname = os.path.join(cif_output_path,
                                  f"{data_prefix}_{atoms_str}_rank{i}_sym{stab_size}.cif")
            smap = {str(j): elem for j, elem in enumerate(res_data['coloring'])}
            try:
                assembler.assemble_and_save(fname, smap, graph_data, forced_edge=forced_edge)
                cif_paths.append(fname)
            except Exception as e:
                print(f"  [Error] Failed to assemble {fname}: {e}")

    out_path = generate_output_filename(data_prefix, all_ligands)
    final_output = {
        "meta": {"topology": data_prefix, "group_size": group_size, "verification": stats},
        "structures": results,
    }
    with open(out_path, 'w') as f:
        json.dump(final_output, f, indent=2)
    print(f"\nResults saved to: {out_path}")

    return cif_paths, out_path, final_output


def main():
    import argparse
    import resource
    import time as time_module

    parser = argparse.ArgumentParser(description="MOF structure enumeration")
    parser.add_argument("--benchmark-json", type=str, default=None,
                        help="If set, write a benchmark JSON record to this path after enumeration")
    parser.add_argument("--disable-budget", action="store_true", default=False,
                        help="Disable all enumeration budget limits for exhaustive runs")
    args = parser.parse_args()

    # ENUM_DISABLE_BUDGET=1 env var also disables budget
    if os.environ.get("ENUM_DISABLE_BUDGET") == "1":
        args.disable_budget = True

    # ================= 配置区 =================
    DATA_PREFIX = "vba"
    ATOM_CHOICES = ["4c_B_Ch", "4c_Cd_1_Ch"]
    # ATOM_CHOICES = ["4c_Cu_1_Ch", "4c_B_Ch", "3c_B_Ch"]
    FORCED_EDGE = "1B_4H_Ch"

    # ── CIF 生成模式 ──────────────────────────────────────────────────────────
    SAMPLE_MODE = True

    MAX_CIFS    = 5
    N_PER_CAT   = 20
    HIGH_RHO    = 0.75
    LOW_RHO     = 0.25
    TASK_TIME_LIMIT = 60
    MAX_TOTAL_RESULTS = 100_000
    # =========================================

    # ── Budget override ──────────────────────────────────────────────────────
    if args.disable_budget:
        enum_task_time_limit = None
        enum_max_total = None
        enum_orbit_limit = float('inf')
        enum_new_limit = None
    else:
        enum_task_time_limit = TASK_TIME_LIMIT
        enum_max_total = MAX_TOTAL_RESULTS
        enum_orbit_limit = None   # None = use module-constant default
        enum_new_limit = None     # None = use module-constant default

    graph_path = os.path.join("data", DATA_PREFIX, f"{DATA_PREFIX}_graph.json")
    if not os.path.exists(graph_path):
        print(f"Error: Graph file not found at {graph_path}")
        return

    print(f"--- Loading Topology: {DATA_PREFIX} ---")

    graph_data = load_graph(graph_path)
    graph_data = enrich_graph_data_with_topology(graph_data)

    # ---- 拓扑需求自动分析 ----
    topology_requirements = calculate_topology_requirements(graph_data)
    print(f"  [Topology] Requirements: {topology_requirements}")

    num_nodes = graph_data["n_nodes"]
    node_degrees = [n.get('degree', 0) for n in graph_data['nodes']]

    # ---- 将 ATOM_CHOICES 转换为 user_selected_ligands ----
    user_selected_ligands = {}
    for atom in ATOM_CHOICES:
        deg = int(atom.split('c')[0])
        user_selected_ligands.setdefault(deg, []).append(atom)
    all_ligands = sorted(set(
        lig for ligs in user_selected_ligands.values() for lig in ligs
    ))
    print(f"  [Ligands] User selected: {user_selected_ligands}")

    # Step 1: 对称性分析
    print("\n--- Step 1: Symmetry Analysis ---")
    sym_engine = SymmetryEngine(graph_data, verbose=True)
    perms = sym_engine.compute_automorphism_group()
    group_size = len(perms)
    print(f"  Automorphism Group Size |G|: {group_size}")

    # Step 2: 结构枚举 (含早期剪枝)
    print("\n--- Step 2: Structure Enumeration ---")
    enumerator = MOFEnumerator(perms, num_nodes, verbose=True,
                               topology_requirements=topology_requirements,
                               user_selected_ligands=user_selected_ligands,
                               node_degrees=node_degrees,
                               task_time_limit=enum_task_time_limit,
                               max_total_results=enum_max_total,
                               identity_orbit_limit=enum_orbit_limit,
                               identity_new_limit=enum_new_limit)

    # ── Benchmark timing/memory (only wraps enumerate()) ──────────────────────
    _bench_t0 = time_module.perf_counter()
    _rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    results = enumerator.enumerate(all_ligands)
    _rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    _bench_t1 = time_module.perf_counter()
    _bench_time_s = _bench_t1 - _bench_t0
    # macOS ru_maxrss is bytes; compute peak as post-call high-watermark
    _bench_peak_mb = _rss_after / 1_048_576

    print(f"  Found {len(results)} unique structures (early-pruned by degree & count constraints).")

    # Step 3: 验证
    print("\n--- Step 3: Verification ---")
    verifier = BridgeVerifier(group_size)
    is_verified, stats = verifier.verify(results, num_nodes, all_ligands, None,
                                          topology_requirements=topology_requirements,
                                          user_selected_ligands=user_selected_ligands)

    # ── Benchmark JSON output (only when --benchmark-json is set) ─────────────
    if args.benchmark_json:
        _ligand_parts = []
        for deg in sorted(user_selected_ligands.keys()):
            _ligand_parts.append(f"C{deg}={len(user_selected_ligands[deg])}")
        _ligand_pool = ",".join(_ligand_parts)
        _bench_record = {
            "topology": DATA_PREFIX,
            "ligand_pool": _ligand_pool,
            "orbits": len(results),
            "time_s": round(_bench_time_s, 6),
            "peak_mem_mb": round(_bench_peak_mb, 2),
            "exhaustive": bool(args.disable_budget),
            "budget_triggered": not is_verified and not args.disable_budget,
        }
        _bench_path = os.path.abspath(args.benchmark_json)
        os.makedirs(os.path.dirname(_bench_path) or ".", exist_ok=True)
        with open(_bench_path, 'w') as _f:
            json.dump(_bench_record, _f, indent=2)
        print(f"  [Benchmark] Wrote {_bench_path}")

    # Step 4: 生成 CIF
    print("\n--- Step 4: Assembler (Generate CIFs) ---")

    atoms_str = "".join(sorted(all_ligands))
    cif_subdir = f"{DATA_PREFIX}_{atoms_str}_cifs"
    cif_output_path = os.path.join("results", cif_subdir)
    os.makedirs(cif_output_path, exist_ok=True)

    assembler = StructureAssembler(output_dir=cif_output_path)

    if SAMPLE_MODE:
        # ── 分层抽样模式 ────────────────────────────────────────────────────
        print(f"  [采样模式] ρ 阈值: High≥{HIGH_RHO} / Med ({LOW_RHO},{HIGH_RHO}) / Low≤{LOW_RHO}")
        print(f"  [采样模式] 每类最多 {N_PER_CAT} 个 CIF，共最多 {N_PER_CAT * 3} 个 …")
        sampled = sample_by_symmetry(results, group_size,
                                      n_per_cat=N_PER_CAT,
                                      high_rho=HIGH_RHO, low_rho=LOW_RHO)
        for i, (res_data, cat, rho) in enumerate(sampled):
            stab = res_data['stabilizer_size']
            fname = os.path.join(
                cif_output_path,
                f"{DATA_PREFIX}_{atoms_str}_{cat}_sym{stab}_rho{rho:.3f}_{i:03d}.cif")
            smap = {str(j): elem for j, elem in enumerate(res_data['coloring'])}
            try:
                assembler.assemble_and_save(fname, smap, graph_data, forced_edge=FORCED_EDGE)
                print(f"    [{cat}] {os.path.basename(fname)}")
            except Exception as e:
                print(f"  [Error] {os.path.basename(fname)}: {e}")

    else:
        # ── 原始 Top-K 模式 ─────────────────────────────────────────────────
        results.sort(key=lambda x: x['stabilizer_size'], reverse=True)
        print(f"  Generating CIFs for Top-{min(len(results), MAX_CIFS)} Highest Symmetry structures...")
        for i, res_data in enumerate(results[:MAX_CIFS]):
            stab_size = res_data['stabilizer_size']
            fname = os.path.join(cif_output_path,
                                  f"{DATA_PREFIX}_{atoms_str}_rank{i}_sym{stab_size}.cif")
            smap = {str(j): elem for j, elem in enumerate(res_data['coloring'])}
            try:
                assembler.assemble_and_save(fname, smap, graph_data, forced_edge=FORCED_EDGE)
            except Exception as e:
                print(f"  [Error] Failed to assemble {fname}: {e}")

    out_path = generate_output_filename(DATA_PREFIX, all_ligands)
    final_output = {
        "meta": {"topology": DATA_PREFIX, "group_size": group_size, "verification": stats},
        "structures": results
    }
    with open(out_path, 'w') as f:
        json.dump(final_output, f, indent=2)
    print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()