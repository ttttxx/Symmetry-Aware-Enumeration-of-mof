
import os, sys, csv, random, re, time, warnings, tempfile, argparse
import multiprocessing as mp
import numpy as np
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

warnings.filterwarnings('ignore')

current_dir  = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.insert(0, project_root)
os.chdir(project_root)

from main import (load_graph, calculate_topology_requirements,
                  enrich_graph_data_with_topology)
from symmetry  import SymmetryEngine
from enumerator import MOFEnumerator
from assembler  import StructureAssembler

results_root = os.path.join(project_root, "results", "structure_selection")
cifs_root    = os.path.join(project_root, "results")

EXPERIMENTS = [
    ("ffc",
     {3: ["3c_B_Ch", "3c_triazine_Ch"],
      4: ["4c_Cu_1_Ch", "4c_B_Ch"]},
     "3c×2 + 4c×2", "enumerate", None),
    ("ffc",
     {3: ["3c_B_Ch", "3c_triazine_Ch", "3c_Cu_imadazoline_Ch"],
      4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"]},
     "3c×3 + 4c×3", "enumerate", None),

    ("fjh",
     {3: ["3c_B_Ch"],
      4: ["4c_Cu_1_Ch", "4c_B_Ch"]},
     "3c×1 + 4c×2", "enumerate", None),
    ("fjh",
    ["fjh_3c_B_Ch3c_triazine_Ch4c_B_Ch4c_Cd_1_Ch4c_Cu_1_Ch_cifs"],
    "3c×2 + 4c×3", "cif_dirs",
     {3: ["3c_B_Ch", "3c_triazine_Ch"],
      4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"]}
     ),

    ("vba",
     ["vba_4c_B_Ch4c_Cu_1_Ch_cifs"],
     "4c×2 (pre-gen)", "cif_dirs",
     {4: ["4c_B_Ch", "4c_Cu_1_Ch"]}),          # fallback

    ("vba",
     ["vba_4c_B_Ch4c_Cd_1_Ch4c_Cu_1_Ch_cifs"],
     "4c×3 (pre-gen)", "cif_dirs",
     {4: ["4c_B_Ch", "4c_Cd_1_Ch", "4c_Cu_1_Ch"]}),                       # fallback
]

ENUM_TIMEOUT  = 600
MAX_HIGH      = 50
MAX_MED       = 12
MAX_LOW       = 12
MAX_ATOMS     = 3000
TORCH_THREADS = 4

E_MIN = -200
E_MAX =  1000.0

VBA_FALLBACK_MAX = 5

MIN_POINTS_ANNOTATE = 3

HIGH_RHO = 0.75   # ρ >= HIGH_RHO  → High
LOW_RHO  = 0.25   # ρ <= LOW_RHO   → Low  |  (LOW_RHO, HIGH_RHO) → Med

STYLE = {
    'High': dict(color='#C62828', marker='o', s=90,  alpha=0.92, zorder=5),
    'Med':  dict(color='#1565C0', marker='^', s=65,  alpha=0.80, zorder=4),
    'Low':  dict(color='#78909C', marker='D', s=45,  alpha=0.65, zorder=3),
}
ZONE_COLOR = {
    'High': '#FFCDD2',  # light red
    'Med':  '#BBDEFB',  # light blue
    'Low':  '#ECEFF1',  # light grey-blue
}

_chgnet_model = None

def get_chgnet():
    global _chgnet_model
    if _chgnet_model is None:
        import torch
        torch.set_num_threads(TORCH_THREADS)
        from chgnet.model import CHGNet
        print(f">>> Loading CHGNet  (PyTorch 线程上限: {TORCH_THREADS}) …")
        _chgnet_model = CHGNet.load()
        print(">>> CHGNet ready.")
    return _chgnet_model

def _enum_worker(queue, topology, user_selected_ligands):
    path = os.path.join(project_root, "data", topology, f"{topology}_graph.json")
    g    = load_graph(path)
    g    = enrich_graph_data_with_topology(g)
    reqs    = calculate_topology_requirements(g)
    n_nodes = g['n_nodes']
    ndeg    = [nd.get('degree', 0) for nd in g['nodes']]
    ligands = sorted({l for ligs in user_selected_ligands.values() for l in ligs})
    perms   = SymmetryEngine(g, verbose=False).compute_automorphism_group()
    enum    = MOFEnumerator(perms, n_nodes, verbose=False,
                            topology_requirements=reqs,
                            user_selected_ligands=user_selected_ligands,
                            node_degrees=ndeg)
    results = enum.enumerate(ligands)
    queue.put((results, len(perms), g))

def _enumerate_with_timeout(topology, user_selected_ligands):
    ctx   = mp.get_context('spawn')
    queue = ctx.Queue()
    p = ctx.Process(target=_enum_worker,
                    args=(queue, topology, user_selected_ligands), daemon=True)
    p.start()
    p.join(timeout=ENUM_TIMEOUT)
    timed_out = p.is_alive()
    if timed_out:
        p.terminate(); p.join()
    if not queue.empty():
        results, group_size, g = queue.get_nowait()
    else:
        results, group_size, g = [], 1, {}
    print(f"  枚举: {len(results)} 个结构  |Γ|={group_size}"
          + (" [超时截断]" if timed_out else ""))
    return results, group_size, g

def _rho_category(rho):
    if rho >= HIGH_RHO: return 'High'
    if rho <= LOW_RHO:  return 'Low'
    return 'Med'

def _stratified_sample(results, group_size, max_per_cat=None):
    pool = {'High': [], 'Med': [], 'Low': []}
    for r in results:
        cat = _rho_category(r['stabilizer_size'] / group_size)
        pool[cat].append(r)

    limits = max_per_cat or {'High': MAX_HIGH, 'Med': MAX_MED, 'Low': MAX_LOW}
    sampled = {
        cat: (lst if len(lst) <= limits[cat] else random.sample(lst, limits[cat]))
        for cat, lst in pool.items()
    }
    for cat, lst in sampled.items():
        print(f"    {cat}: pool={len(pool[cat])}  sampled={len(lst)}")
    return sampled

def _assemble_and_evaluate(coloring, graph_data, assembler, tag):
    from pymatgen.core import Structure
    fd, tmp = tempfile.mkstemp(suffix='.cif', prefix=f'chgnet_{tag}_')
    os.close(fd)
    try:
        smap = {str(i): atom for i, atom in enumerate(coloring)}
        assembler.assemble_and_save(tmp, smap, graph_data)
        struct = Structure.from_file(tmp)
        if len(struct) > MAX_ATOMS:
            print(f"      [熔断-原子数] {len(struct)} > {MAX_ATOMS}")
            return None
        E = float(get_chgnet().predict_structure(struct)['e'])
        if not (E_MIN <= E <= E_MAX):
            print(f"      [熔断-能量异常] E={E:.2f} eV/atom 超出 [{E_MIN}, {E_MAX}]，"
                  f"结构很可能由 CIF 组装错误导致，已丢弃")
            return None
        return E
    except Exception as ex:
        print(f"      [Skip] {ex}")
        return None
    finally:
        for f in [tmp, tmp.replace('.cif', '.log')]:
            try:
                if os.path.exists(f): os.remove(f)
            except OSError:
                pass

def run_enumerate_experiment(exp_idx, topology, user_selected_ligands, label,
                              max_per_cat=None):
    print(f"\n{'='*62}")
    print(f"[实验 {exp_idx+1}]  {topology.upper()} | {label}  [枚举模式]")
    print('='*62)
    t0 = time.time()

    results, group_size, graph_data = _enumerate_with_timeout(
        topology, user_selected_ligands)
    if not results:
        return [], _empty_meta(label, topology, group_size)

    sampled = _stratified_sample(results, group_size, max_per_cat=max_per_cat)
    total   = sum(len(v) for v in sampled.values())
    print(f"  → 需评估 {total} 个结构")

    assembler   = StructureAssembler(output_dir=tempfile.gettempdir())
    data_points = []
    idx         = 0

    for cat, structs in [('High', sampled['High']),
                          ('Med',  sampled['Med']),
                          ('Low',  sampled['Low'])]:
        print(f"\n  >> {cat} ({len(structs)} 个)")
        for r in structs:
            stab = r['stabilizer_size']
            rho  = stab / group_size
            tc   = time.time()
            E    = _assemble_and_evaluate(r['coloring'], graph_data, assembler,
                                          f"{topology}_{exp_idx}_{idx}")
            idx += 1
            if E is not None:
                data_points.append(dict(rho=rho, stab=stab, energy=E, category=cat))
                print(f"    [{idx:3d}] stab={stab:3d} ρ={rho:.4f}  "
                      f"E={E:.4f} eV/atom  [{time.time()-tc:.0f}s]")
            else:
                print(f"    [{idx:3d}] stab={stab:3d} ρ={rho:.4f}  "
                      f"SKIP [{time.time()-tc:.0f}s]")

    meta = dict(label=label, topology=topology, group_size=group_size)
    print(f"\n  完成: {len(data_points)}/{total} 点  "
          f"耗时 {(time.time()-t0)/60:.1f} min")
    return data_points, meta

def run_cif_dirs_experiment(exp_idx, topology, cif_dir_names, label,
                             fallback_ligands=None, vba_fallback_max=None):
    print(f"\n{'='*62}")
    print(f"[实验 {exp_idx+1}]  {topology.upper()} | {label}  [预生成CIF模式]")
    print('='*62)
    from pymatgen.core import Structure

    cif_entries = []
    for dname in cif_dir_names:
        d = os.path.join(cifs_root, dname)
        if not os.path.exists(d):
            print(f"  [Warning] CIF目录不存在: {d}")
            continue
        for f in sorted(os.listdir(d)):
            if not f.endswith('.cif'): continue
            m = re.search(r'_sym(\d+)', f)
            sym = int(m.group(1)) if m else 1
            cif_entries.append((os.path.join(d, f), sym))

    if not cif_entries:
        if fallback_ligands is not None:
            n = vba_fallback_max or VBA_FALLBACK_MAX
            est_h = n * 3 * 800 / 3600
            print(f"\n  [Fallback] 所有预生成目录均不存在。")
            print(f"  [Fallback] 切换枚举模式，每类最多 {n} 个结构。")
            print(f"  [Fallback] 注意：{topology.upper()} CIF 组装约 800 s/个，")
            print(f"             预计总耗时 {n}×3×800s ≈ {est_h:.1f} 小时。")
            print(f"             可用 --vba-n 调整每类数量。")
            return run_enumerate_experiment(
                exp_idx, topology, fallback_ligands,
                label + " [fallback]",
                max_per_cat={'High': n, 'Med': n, 'Low': n})
        else:
            print("  [Error] 无CIF文件且未配置 fallback_ligands，跳过此实验。")
            return [], _empty_meta(label, topology, 16)

    group_size  = 16
    data_points = []
    model       = get_chgnet()

    for i, (cif_path, sym) in enumerate(cif_entries):
        rho   = sym / group_size
        cat   = _rho_category(rho)
        fname = os.path.basename(cif_path)
        print(f"  [{i+1:2d}] {fname}  sym={sym} ρ={rho:.4f} cat={cat} …",
              end='', flush=True)
        try:
            struct = Structure.from_file(cif_path)
            if len(struct) > MAX_ATOMS:
                print(f"  [熔断-原子数] {len(struct)}")
                continue
            E = float(model.predict_structure(struct)['e'])
            if not (E_MIN <= E <= E_MAX):
                print(f"  [熔断-能量异常] E={E:.2f}")
                continue
            data_points.append(dict(rho=rho, stab=sym, energy=E, category=cat))
            print(f"  E={E:.4f} eV/atom")
        except Exception as ex:
            print(f"  FAIL: {ex}")

    meta = dict(label=label, topology=topology, group_size=group_size)
    print(f"\n  完成: {len(data_points)} 点")
    return data_points, meta

def _empty_meta(label, topology, group_size):
    return dict(label=label, topology=topology, group_size=group_size)

def save_csv(data_points, meta, filepath):
    with open(filepath, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["Filename", "Label", "Symmetry_Order",
                    "Energy_eV_per_atom", "Symmetry_Index_Rho", "Type"])
        for i, p in enumerate(data_points):
            w.writerow([
                f"{meta['topology']}_struct_{i:04d}.cif",
                f"{meta['topology']}_{meta['label']}",
                p['stab'],
                f"{p['energy']:.6f}",
                f"{p['rho']:.6f}",
                p['category'],
            ])
    print(f"  [CSV] {filepath}")

def _setup_style():
    """Apply Nature / ACS journal–style rcParams globally."""
    plt.rcParams.update({
        'font.family':           'sans-serif',
        'font.sans-serif':       ['Arial', 'DejaVu Sans', 'Helvetica'],
        'font.size':             9,
        'axes.linewidth':        0.9,
        'axes.labelsize':        9,
        'axes.titlesize':        10,
        'axes.titleweight':      'bold',
        'axes.spines.top':       False,
        'axes.spines.right':     False,
        'xtick.major.width':     0.9,
        'ytick.major.width':     0.9,
        'xtick.minor.width':     0.6,
        'ytick.minor.width':     0.6,
        'xtick.direction':       'in',
        'ytick.direction':       'in',
        'xtick.major.size':      4,
        'ytick.major.size':      4,
        'xtick.minor.size':      2.5,
        'ytick.minor.size':      2.5,
        'legend.fontsize':       7.5,
        'legend.framealpha':     0.92,
        'legend.edgecolor':      '#cccccc',
        'legend.handlelength':   1.5,
        'figure.dpi':            150,
        'savefig.dpi':           300,
        'savefig.bbox':          'tight',
    })

def _fmt(val, _pos=None):
    return f"{val:.4f}".rstrip('0').rstrip('.')

def _draw_zone_bands(ax, data_points):
    rhos = [p['rho'] for p in data_points]
    x_lo = min(rhos) * 0.7
    x_hi = max(rhos) * 1.4

    bands = [
        (max(x_lo, HIGH_RHO), x_hi,          'High'),
        (max(x_lo, LOW_RHO),  min(x_hi, HIGH_RHO), 'Med'),
        (x_lo,                 min(x_hi, LOW_RHO),  'Low'),
    ]
    for x0, x1, cat in bands:
        if x0 < x1:
            ax.axvspan(x0, x1, alpha=0.07,
                       color=ZONE_COLOR[cat], zorder=0, lw=0)

    for thresh in [HIGH_RHO, LOW_RHO]:
        if min(rhos) < thresh < max(rhos):
            ax.axvline(thresh, color='#bbbbbb', lw=0.7, ls=':', zorder=1)

def _draw_trend(ax, data_points):
    if len(data_points) < 4:
        return
    log_rho  = np.array([np.log10(p['rho'])    for p in data_points])
    energies = np.array([p['energy']            for p in data_points])
    slope, intercept, r_val, p_val, _ = stats.linregress(log_rho, energies)

    x_trend = np.logspace(log_rho.min(), log_rho.max(), 300)
    y_trend = slope * np.log10(x_trend) + intercept
    ax.plot(x_trend, y_trend, '--', color='#444444',
            lw=1.1, alpha=0.50, zorder=2, label=f'Linear fit')

    sig = '***' if p_val < 0.001 else ('**' if p_val < 0.01 else
          ('*' if p_val < 0.05 else 'n.s.'))
    ax.text(0.97, 0.96,
            f'$r$ = {r_val:.2f}{sig}',
            transform=ax.transAxes,
            ha='right', va='top',
            fontsize=8, style='italic',
            bbox=dict(boxstyle='round,pad=0.25',
                      facecolor='white', edgecolor='#cccccc',
                      alpha=0.85, lw=0.6))

def _annotate_ground_state(ax, data_points):
    best  = min(data_points, key=lambda p: p['energy'])
    enes  = [p['energy'] for p in data_points]
    y_range = max(enes) - min(enes)

    y_text = best['energy'] + y_range * 0.08

    ax.annotate(
        f'Min $E$\n($\\rho$ = {_fmt(best["rho"])})',
        xy=(best['rho'], best['energy']),
        xytext=(best['rho'], y_text),
        fontsize=8,
        ha='center',
        va='bottom',
        fontweight='bold',
        color='#111111',
        arrowprops=dict(
            arrowstyle='-|>',
            color='#111111',
            lw=1.5,
            mutation_scale=10,
        ),
        zorder=10,
    )

def _draw_ax(ax, data_points, meta):
    if not data_points:
        ax.text(0.5, 0.5, "No data\n(CIF generation\nnot feasible)",
                ha='center', va='center', transform=ax.transAxes,
                fontsize=10, color='#888888')
        ax.set_title(f"{meta['topology'].upper()} | {meta['label']}")
        return

    _draw_zone_bands(ax, data_points)

    groups = {'High': [], 'Med': [], 'Low': []}
    for p in data_points:
        groups[p['category']].append(p)

    legend_handles = []
    for cat in ['Low', 'Med', 'High']:
        pts = groups[cat]
        if not pts: continue
        x  = [p['rho']    for p in pts]
        y  = [p['energy'] for p in pts]
        st = STYLE[cat]
        ax.scatter(x, y,
                   c=st['color'], marker=st['marker'],
                   s=st['s'],     alpha=st['alpha'],
                   edgecolors='white', linewidths=0.6,
                   zorder=st['zorder'])

        if cat == 'High':
            label_txt = fr'High symmetry  ($\rho \geq {HIGH_RHO}$)'
        elif cat == 'Med':
            label_txt = fr'Med. symmetry  (${LOW_RHO} < \rho < {HIGH_RHO}$)'
        else:
            label_txt = fr'Low symmetry  ($\rho \leq {LOW_RHO}$)'

        legend_handles.append(
            Line2D([0], [0],
                   marker=st['marker'], color='w',
                   markerfacecolor=st['color'],
                   markeredgecolor='white',
                   markersize=7 if cat == 'High' else (6 if cat == 'Med' else 5.5),
                   label=label_txt))

    if len(data_points) >= MIN_POINTS_ANNOTATE:
        _annotate_ground_state(ax, data_points)

    ax.set_title(f"{meta['topology'].upper()}  |  {meta['label']}")
    ax.set_xlabel(r"Symmetry Index  $\rho(c) = |\mathrm{Stab}(c)|\,/\,|\Gamma|$",
                  fontsize=8.5)
    ax.set_ylabel("CHGNet Energy  (eV atom$^{-1}$)", fontsize=8.5)
    ax.set_xscale('log')

    ux = sorted({p['rho'] for p in data_points})
    ax.set_xticks(ux)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(_fmt))
    ax.tick_params(axis='x', labelsize=7.5, rotation=20)
    ax.tick_params(axis='y', labelsize=8)

    ax.yaxis.grid(True, which='major', ls='--', lw=0.5, alpha=0.35, color='#aaaaaa')
    ax.set_axisbelow(True)

    ax.legend(handles=legend_handles,
              loc='best', framealpha=0.92,
              fontsize=7.2, borderpad=0.6, labelspacing=0.45)

def _safe_label(label):
    return (label.replace(' ', '').replace(',', '_')
                 .replace('×', 'x').replace('+', '')
                 .replace('(', '').replace(')', ''))

def plot_single(data_points, meta, exp_idx):
    _setup_style()
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 4.8))
    fig.suptitle(
        "CHGNet Symmetry–Stability Correlation\n"
        f"{meta['topology'].upper()}  |  {meta['label']}",
        fontsize=11, fontweight='bold', y=1.01,
    )
    _draw_ax(ax, data_points, meta)
    plt.tight_layout()

    safe = _safe_label(meta['label'])
    out  = os.path.join(results_root,
                        f"{meta['topology']}_{safe}_exp{exp_idx+1}.png")
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [单图] {out}")

def plot_combined(all_data):
    _setup_style()
    fig, axes = plt.subplots(3, 2, figsize=(13, 14))
    fig.suptitle(
        "CHGNet Symmetry–Stability Correlation\n"
        r"$\rho(c) = |\mathrm{Stab}(c)|\,/\,|\Gamma|$   vs   "
        r"Energy (eV atom$^{-1}$)",
        fontsize=13, fontweight='bold', y=0.998,
    )
    for (data_points, meta), ax in zip(all_data, axes.flat):
        _draw_ax(ax, data_points, meta)

    plt.subplots_adjust(hspace=0.52, wspace=0.32)
    out = os.path.join(results_root, "chgnet_combined_analysis.png")
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n[总图] {out}")

def main():
    parser = argparse.ArgumentParser(
        description="CHGNet MOF symmetry-stability analysis")
    parser.add_argument(
        '--exp', type=int, default=None, metavar='N',
        help='只跑第 N 个实验（0-based，共 0-5）。不指定则跑全部 6 个。')
    parser.add_argument(
        '--vba-n', type=int, default=None, metavar='N',
        help=(f'vba fallback 时每类最多评估 N 个结构（默认 {VBA_FALLBACK_MAX}）。'
              f'注意每个 CIF ~800 s，3×N 个结构约需 N×40 分钟。'))
    args = parser.parse_args()

    if args.exp is not None and not (0 <= args.exp < len(EXPERIMENTS)):
        parser.error(f"--exp 必须在 0 到 {len(EXPERIMENTS)-1} 之间")

    os.makedirs(results_root, exist_ok=True)
    random.seed(42)
    get_chgnet()

    exp_list = (list(enumerate(EXPERIMENTS)) if args.exp is None
                else [(args.exp, EXPERIMENTS[args.exp])])

    all_data = []
    t_start  = time.time()

    for idx, exp_entry in exp_list:
        topo, config, label, mode = exp_entry[:4]
        fallback_ligands = exp_entry[4] if len(exp_entry) > 4 else None

        if mode == 'enumerate':
            dp, meta = run_enumerate_experiment(idx, topo, config, label)
        else:
            dp, meta = run_cif_dirs_experiment(
                idx, topo, config, label,
                fallback_ligands=fallback_ligands,
                vba_fallback_max=args.vba_n)

        all_data.append((dp, meta))
        plot_single(dp, meta, idx)

        if dp:
            safe  = _safe_label(label)
            csv_p = os.path.join(results_root,
                                 f"{topo}_{safe}_chgnet_data.csv")
            save_csv(dp, meta, csv_p)

    if args.exp is None:
        plot_combined(all_data)

    elapsed = (time.time() - t_start) / 60
    print(f"\n>>> 完成。总耗时: {elapsed:.1f} 分钟")

if __name__ == "__main__":
    mp.freeze_support()
    main()