# !/usr/bin/env python3

import os, sys, time, random, pickle
from typing import List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from main import load_graph, calculate_topology_requirements
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator

sns.set_style("whitegrid", {
    "ytick.major.size": 0,
    "xtick.major.size": 0,
    "grid.color": "#f0f0f0",
    "grid.linestyle": "-",
})
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.8,
    "lines.linewidth": 2.2,
    "legend.fontsize": 9.5,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

PALETTE = sns.color_palette("colorblind", 4)
COLOR_OURS      = PALETTE[0]    # stratified enumeration (blue)
COLOR_BASELINE  = PALETTE[3]    # random baseline (reddish)
COLOR_OURS_BS   = "#6baed6"     # lighter blue for Burnside
COLOR_BASE_BS   = "#fd8d6c"     # lighter red for baseline Burnside
COLOR_CUTOFF    = "#444444"
COLOR_HS_LINE   = "#2ca02c"

class TrackedEnumerator(MOFEnumerator):
    """Thin subclass that appends ``stabilizer_size`` to ``emission_log``
    inside ``_try_add_result``, preserving the *discovery order* that the
    original ``enumerate()`` discards when it sorts by |Stab| at the end."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.emission_log: List[int] = []

    def _try_add_result(self, coloring):
        raw_tuple = tuple(coloring)
        canonical_tup = self._get_canonical(raw_tuple)
        if self._trie_lookup(canonical_tup) is not None:
            return
        stab_size, _ = self._compute_stabilizer_and_orbit(canonical_tup)
        cid = len(self.classes)
        self._trie_insert(canonical_tup, cid, stab_size)
        self.classes[cid] = {
            "id": cid,
            "coloring": list(canonical_tup),
            "stabilizer_size": stab_size,
            "orbit_size": len(self.perms) // stab_size,
        }
        self.emission_log.append(stab_size)

def _theoretical_orbits(topo_req: dict, ligands: dict,
                        fallback_atoms: list, g_size: int) -> float:
    """|X| / |Γ| = Π_d |L_d|^{count_d} / |G|"""
    total = 1
    for deg, cnt in topo_req.items():
        n = len(ligands.get(deg, fallback_atoms))
        total *= (n ** cnt)
    return total / max(1, g_size)

def _ligand_key(ligands: dict) -> str:
    """Short hash string for the ligand dict."""
    parts = []
    for deg in sorted(ligands):
        atoms = "-".join(sorted(ligands[deg]))
        parts.append(f"d{deg}_{atoms}")
    return "__".join(parts)

CACHE_DIR = os.path.join(PROJECT_ROOT, "results", "cutoff", "_cache")

def _cache(topology: str, ligands: dict, suffix: str) -> str:
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR,
                        f"{topology}__{_ligand_key(ligands)}__{suffix}.pkl")

TASK_TIME_LIMIT = 60            # seconds per subgroup task
MAX_TOTAL_RESULTS = 500_000     # global structure cap

def run_stratified(topology: str,
                   ligands: dict,
                   *,
                   task_time_limit: float = TASK_TIME_LIMIT,
                   max_total_results: int = MAX_TOTAL_RESULTS,
                   force_rerun: bool = False,
                   verbose: bool = True) -> dict:
    """Stratified enumeration → cumulative curves (cached to disk).

    Each subgroup task is capped at ``task_time_limit`` seconds (60 by
    default).  Enumeration halts globally once ``max_total_results``
    structures are found.  Both guards prevent silent hangs on large
    configuration spaces.
    """
    cpath = _cache(topology, ligands, "stratified")
    if not force_rerun and os.path.exists(cpath):
        if verbose:
            print(f"  [Cache] Load stratified ← {os.path.basename(cpath)}")
        with open(cpath, "rb") as f:
            return pickle.load(f)

    graph_path = os.path.join(PROJECT_ROOT, "data", topology,
                              f"{topology}_graph.json")
    graph_data = load_graph(graph_path)
    all_ligs = sorted({l for ligs in ligands.values() for l in ligs})
    node_degrees = [n.get("degree", 0) for n in graph_data["nodes"]]
    topo_req = calculate_topology_requirements(graph_data)

    sym = SymmetryEngine(graph_data, verbose=False)
    perms = sym.compute_automorphism_group()
    G = len(perms)

    if verbose:
        print(f"  [Safety] task_time_limit={task_time_limit}s / task, "
              f"max_total_results={max_total_results}")

    t0 = time.monotonic()
    enum = TrackedEnumerator(
        perms, graph_data["n_nodes"],
        verbose=verbose,
        topology_requirements=topo_req,
        user_selected_ligands=ligands,
        node_degrees=node_degrees,
        task_time_limit=task_time_limit,        # ← 60 s per subgroup
        max_total_results=max_total_results,    # ← 500 k global cap
    )
    results = enum.enumerate(all_ligs)
    wall = time.monotonic() - t0

    theoretical = _theoretical_orbits(topo_req, ligands, all_ligs, G)
    rho_all = [r["stabilizer_size"] / G for r in results]
    total_hs = sum(1 for rho in rho_all if rho >= 0.75)

    steps, hs_cap, bs_sum, bs_cov = [], [], [], []
    cum_hs, cum_bs = 0, 0.0
    hs_plateau_step = None

    for idx, stab in enumerate(enum.emission_log):
        rho = stab / G
        cum_hs += (rho >= 0.75)
        cum_bs += 1.0 / stab
        s = idx + 1
        steps.append(s)
        hs_cap.append(cum_hs / total_hs * 100.0 if total_hs else 0.0)
        bs_sum.append(cum_bs)
        bs_cov.append(cum_bs / theoretical * 100.0 if theoretical else 0.0)
        if hs_plateau_step is None and cum_hs >= total_hs:
            hs_plateau_step = s

    if verbose:
        print(f"  [Stratified] {len(results)} orbits in {wall:.1f}s | "
              f"G={G} | HS={total_hs} (plateau @ step {hs_plateau_step}) | "
              f"|X|/|G|={theoretical:.2e} | "
              f"Σ1/|Stab|={cum_bs:.2f}")

    data = {
        "steps": steps, "high_sym_cap_pct": hs_cap,
        "burnside_sum": bs_sum, "burnside_cov_pct": bs_cov,
        "group_size": G, "total_high_sym": total_hs,
        "theoretical_orbits": theoretical, "total_found": len(results),
        "hs_plateau_step": hs_plateau_step, "wall_seconds": wall,
    }
    with open(cpath, "wb") as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    return data

def run_baseline(topology: str,
                 ligands: dict,
                 *,
                 num_samples: int = 150_000,
                 ground_truth: dict,
                 seed: int = 42,
                 force_rerun: bool = False,
                 verbose: bool = True) -> dict:
    """Random-sampling baseline — per‑emission cumulative curves (cached)."""
    cpath = _cache(topology, ligands, "baseline")
    if not force_rerun and os.path.exists(cpath):
        if verbose:
            print(f"  [Cache] Load baseline ← {os.path.basename(cpath)}")
        with open(cpath, "rb") as f:
            return pickle.load(f)

    graph_path = os.path.join(PROJECT_ROOT, "data", topology,
                              f"{topology}_graph.json")
    graph_data = load_graph(graph_path)
    n_nodes = graph_data["n_nodes"]
    sym = SymmetryEngine(graph_data, verbose=False)
    perms = tuple(sym.compute_automorphism_group())
    G = len(perms)

    total_hs = ground_truth["total_high_sym"]
    theoretical = ground_truth["theoretical_orbits"]

    rng = random.Random(seed)
    all_ligs = sorted({l for ligs in ligands.values() for l in ligs})

    def canonicalize(col):
        best = col
        for p in perms:
            cand = tuple(col[i] for i in p)
            if cand < best:
                best = cand
        return best

    def stabilizer(col):
        cnt = 0
        for p in perms:
            if all(col[i] == col[p[i]] for i in range(n_nodes)):
                cnt += 1
        return cnt

    found = set()
    cum_hs, cum_bs = 0, 0.0
    steps, hs_cap, bs_sum, bs_cov = [], [], [], []
    t0 = time.monotonic()

    for k in range(1, num_samples + 1):
        coloring = tuple(rng.choices(all_ligs, k=n_nodes))
        canon = canonicalize(coloring)
        if canon in found:
            continue

        found.add(canon)
        stab = stabilizer(canon)
        rho = stab / G
        cum_hs += (rho >= 0.75)
        cum_bs += 1.0 / stab

        steps.append(len(found))
        hs_cap.append(cum_hs / total_hs * 100.0 if total_hs else 0.0)
        bs_sum.append(cum_bs)
        bs_cov.append(cum_bs / theoretical * 100.0 if theoretical else 0.0)

        if verbose and k % 50_000 == 0:
            elapsed = time.monotonic() - t0
            print(f"    [{k:>7d}/{num_samples}] unique={len(found):>6d}  "
                  f"hs={cum_hs}/{total_hs}  Σ1/|Stab|={cum_bs:.1f}  "
                  f"[{elapsed:.0f}s]")

    wall = time.monotonic() - t0
    if verbose:
        n = len(found)
        print(f"  [Baseline] {num_samples} samples → {n} unique "
              f"({n / num_samples * 100:.2f}%) in {wall:.1f}s | "
              f"hs={cum_hs}/{total_hs}")

    data = {
        "steps": steps, "high_sym_cap_pct": hs_cap,
        "burnside_sum": bs_sum, "burnside_cov_pct": bs_cov,
        "total_found": len(found), "wall_seconds": wall,
    }
    with open(cpath, "wb") as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    return data

def make_figure(strat: dict,
                base: dict,
                topology: str,
                ligand_label: str,
                *,
                budget_cutoff: Optional[int] = None,
                output_dir: Optional[str] = None) -> str:
    """Draw a single dual‑axis panel and save as 600 dpi PDF.

    All text elements (legend, info box, annotations) use relative
    positioning so they adapt to the data range without manual tuning.
    """
    if output_dir is None:
        output_dir = os.path.join(PROJECT_ROOT, "results", "cutoff")
    os.makedirs(output_dir, exist_ok=True)

    G = strat["group_size"]
    theoretical = strat["theoretical_orbits"]

    if budget_cutoff is None:
        budget_cutoff = strat.get("hs_plateau_step", len(strat["steps"]))

    x_max = max(max(strat["steps"]), max(base["steps"]))

    fig, ax1 = plt.subplots(figsize=(8.2, 5.8))

    # Push grid lines to the very back so they never occlude legend / text
    ax1.set_axisbelow(True)
    ax1.grid(zorder=0.1)

    ax1.plot(strat["steps"], strat["high_sym_cap_pct"],
             color=COLOR_OURS, linewidth=2.5, linestyle="-", alpha=0.95,
             zorder=4,
             label="Stratified Canonical Enum. (ours) — HS capture")
    ax1.plot(base["steps"], base["high_sym_cap_pct"],
             color=COLOR_BASELINE, linewidth=2.0, linestyle="--", alpha=0.90,
             zorder=3,
             label="Unordered DFS (baseline) — HS capture")

    ax1.set_xscale("log")
    ax1.set_xlabel("Emitted Orbits Index (log scale)", fontweight="bold")
    ax1.set_ylabel("High-Symmetry Capture Rate (%)", fontweight="bold",
                   color=COLOR_OURS)
    ax1.tick_params(axis="y", labelcolor=COLOR_OURS)
    ax1.set_ylim(-5, 112)
    ax1.yaxis.set_major_formatter(mtick.PercentFormatter())

    ax2 = ax1.twinx()
    ax2.set_axisbelow(True)
    ax2.grid(zorder=0.1)
    ax2.plot(strat["steps"], strat["burnside_sum"],
             color=COLOR_OURS_BS, linewidth=2.2, linestyle="-", alpha=0.70,
             zorder=2,
             label="Stratified — Burnside sum")
    ax2.plot(base["steps"], base["burnside_sum"],
             color=COLOR_BASE_BS, linewidth=2.0, linestyle="--", alpha=0.70,
             zorder=1,
             label="Baseline — Burnside sum")

    ax2.set_ylabel(r"Cumulative Burnside Sum  $\Sigma\,1/|{\rm Stab}|$",
                   fontweight="bold", color="#555555")
    ax2.tick_params(axis="y", labelcolor="#555555")

    ax1.axhline(100.0, color=COLOR_HS_LINE, linestyle=":", linewidth=1.6,
                alpha=0.85, zorder=4.5)
    ax1.annotate("Complete Symmetric\nSubspace",
                 xy=(1.6, 103), fontsize=8.2, color=COLOR_HS_LINE,
                 va="bottom", ha="left", style="italic",
                 annotation_clip=False, zorder=12)

    ax1.axvline(budget_cutoff, color=COLOR_CUTOFF, linestyle="--",
                linewidth=1.8, alpha=0.80, zorder=4.5)

    ax1.annotate(f"Budget Cutoff  ({budget_cutoff} orbits)",
                 xy=(budget_cutoff, 97),
                 xytext=(budget_cutoff * 1.8, 104),
                 fontsize=8.8, fontweight="bold", color=COLOR_CUTOFF,
                 arrowprops=dict(arrowstyle="->", color=COLOR_CUTOFF,
                                 lw=1.2, connectionstyle="arc3,rad=0.25"),
                 annotation_clip=False, zorder=12)

    info_str = (f"|G| = {G}     "
                f"|X|/|G| = {theoretical:.2e}     "
                f"HS targets = {strat['total_high_sym']}     "
                f"HS plateau @ step {strat['hs_plateau_step']}")
    fig.text(0.5, 0.01, info_str,
             fontsize=7.2, ha="center", va="bottom",
             bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                       alpha=0.95, edgecolor="#cccccc", linewidth=0.6),
             zorder=15)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    legend = ax1.legend(lines1 + lines2, labels1 + labels2,
                        loc="best", frameon=True,
                        edgecolor="#999999", facecolor="white",
                        framealpha=1.0,
                        fontsize=8.2, ncol=1,
                        labelspacing=0.4,
                        handlelength=1.9, handletextpad=0.7,
                        borderpad=0.5)
    legend.set_zorder(20)

    ax1.set_title(
        f"Cumulative Discovery under Fixed Budget\n"
        f"{topology}  —  {ligand_label}",
        fontsize=12, fontweight="bold", pad=16, linespacing=1.35)

    ax1.set_xlim(0.8, x_max * 1.08)

    fig.tight_layout(pad=1.6, rect=(0.02, 0.06, 0.98, 0.96))

    safe = (ligand_label.replace(" ", "_").replace(",", "")
            .replace("{", "").replace("}", "").replace("|", "_"))
    fname = f"cutoff_{topology}_{safe}.pdf"
    path = os.path.join(output_dir, fname)
    fig.savefig(path, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"  → Saved {path}")
    return path

def main():
    configs = [
        ("vba", {4: ["4c_Cu_1_Ch", "4c_B_Ch"]}),
        ("vba", {4: ["4c_Cu_1_Ch", "4c_B_Ch", "4c_Cd_1_Ch"]}),
    ]
    baseline_samples = 150_000

    for topology, ligands in configs:
        lig_label = "  ".join(
            f"d{deg}:{{{','.join(atoms)}}}"
            for deg, atoms in sorted(ligands.items())
        )
        print(f"\n{'='*70}")
        print(f"  {topology}  |  {lig_label}")
        print(f"{'='*70}")

        print("\n[1/2] Stratified Canonical Enumeration …")
        strat = run_stratified(topology, ligands, verbose=True)

        print(f"\n[2/2] Baseline Unordered DFS ({baseline_samples} samples) …")
        base = run_baseline(topology, ligands,
                            num_samples=baseline_samples,
                            ground_truth=strat, verbose=True)

        cutoff = strat["hs_plateau_step"] or len(strat["steps"])
        print(f"\n  → Budget cutoff = {cutoff} orbits "
              f"(HS plateau: {strat['total_high_sym']} high-sym targets)")

        print("\n[Plot] Generating figure …")
        make_figure(strat, base, topology, lig_label, budget_cutoff=cutoff)

    print(f"\n{'='*70}")
    print("  All done — figures saved under results/cutoff/")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
