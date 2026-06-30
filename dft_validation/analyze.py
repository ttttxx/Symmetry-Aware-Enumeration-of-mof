import os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from dft_validation import config


RHO_CMAP = "viridis"


def _spearman_for_comp(records: list) -> tuple:
    """
    Spearman rs for (rho, dft_e_per_atom) within a composition.
    Returns (rs, p_val, lowest_at_max_rho: bool, n_converged).
    """
    conv = [r for r in records if r.get("dft_e_per_atom") is not None]
    if len(conv) < 2:
        return (float("nan"), float("nan"), False, len(conv))
    rhos     = [r["rho"] for r in conv]
    energies = [r["dft_e_per_atom"] for r in conv]
    # If rho has no variance Spearman is undefined → report NaN cleanly.
    if len(set(rhos)) < 2:
        rs, p = float("nan"), float("nan")
    else:
        rs, p = spearmanr(rhos, energies)
    min_rec = min(conv, key=lambda r: r["dft_e_per_atom"])
    max_rho = max(rhos)
    return (rs, p, min_rec.get("rho") == max_rho, len(conv))


def get_fjh_gate(comp_results: dict) -> tuple:
    """
    Qualitative gate (kept for backward compat / informational use):
    every composition must have lowest-E at the maximum sampled rho AND rs < 0.
    """
    if not comp_results:
        return None, "no compositions analysed"
    for label, stats in comp_results.items():
        if stats.get("n", 0) < 2:
            return None, f"{label}: insufficient converged DFT points"
        if not stats.get("lowest_at_max_rho", False):
            return False, f"{label}: lowest-energy structure is NOT at max rho"
        rs = stats.get("spearman_rs", 0)
        if rs == rs and rs >= 0:  # not NaN and non-negative
            return False, (f"{label}: Spearman rs={rs:.3f} >= 0 "
                           f"(energy does not decrease with symmetry)")
    return True, "all compositions pass the qualitative criterion"


def analyze_topology(records: list, results_dir: str, topology: str) -> dict:
    """
    Per-composition Spearman for one topology + one rho–E figure.
    Returns {composition_label: {spearman_rs, p_val, lowest_at_max_rho, n}}.
    """
    labels = sorted({r["composition_label"] for r in records
                     if r.get("topology") == topology and not r.get("probe_cu")})
    if not labels:
        print(f"[analyze:{topology}] no compositions to plot")
        return {}

    comp_results = {}
    fig, axes = plt.subplots(1, len(labels), figsize=(6 * len(labels), 5), sharey=False)
    if len(labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, labels):
        subset = [r for r in records
                  if r.get("topology") == topology
                  and r.get("composition_label") == label
                  and r.get("dft_e_per_atom") is not None]
        rs, p, lowest_at_max_rho, n = _spearman_for_comp(subset)
        comp_results[label] = {"spearman_rs": rs, "p_val": p,
                               "lowest_at_max_rho": lowest_at_max_rho, "n": n}

        if subset:
            ax.scatter([r["rho"] for r in subset],
                       [r["dft_e_per_atom"] for r in subset],
                       c=[r["rho"] for r in subset],
                       cmap=RHO_CMAP, s=65, alpha=0.9,
                       edgecolor="black", linewidth=0.4, zorder=3)

        ax.set_xlabel(r"Symmetry Index $\rho$", fontsize=12)
        ax.set_ylabel("DFT Energy (eV/atom)", fontsize=12)
        rs_str = f"{rs:.3f}" if rs == rs else "n/a"
        p_str  = f"{p:.3f}" if p == p else "n/a"
        ax.set_title(f"{topology} {label}  (n={n}, $\\rho_s$={rs_str}, p={p_str})", fontsize=11)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(results_dir, f"rho_energy_{topology}.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[analyze:{topology}] Figure saved to {fig_path}")

    if topology == "ffc":
        _plot_pooled_topology(records, results_dir, topology)
    return comp_results


def _plot_pooled_topology(records: list, results_dir: str, topology: str) -> None:
    """Descriptive pooled plot; not used as the controlled statistical test."""
    pts = [r for r in records
           if r.get("topology") == topology
           and not r.get("probe_cu")
           and r.get("dft_e_per_atom") is not None]
    if not pts:
        return
    rs, p = spearmanr([r["rho"] for r in pts], [r["dft_e_per_atom"] for r in pts])
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter([r["rho"] for r in pts],
               [r["dft_e_per_atom"] for r in pts],
               c=[r["rho"] for r in pts], cmap=RHO_CMAP,
               s=65, alpha=0.9, edgecolor="black", linewidth=0.4)
    ax.set_xlabel(r"Symmetry Index $\rho$", fontsize=12)
    ax.set_ylabel("DFT Energy (eV/atom)", fontsize=12)
    ax.set_title(f"{topology} pooled descriptive plot  (n={len(pts)}, "
                 f"$\\rho_s$={rs:.3f}, p={p:.3f})", fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig_path = os.path.join(results_dir, f"rho_energy_{topology}_pooled.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[analyze:{topology}] Pooled descriptive figure saved to {fig_path}")


def analyze_all(records: list, results_dir: str) -> dict:
    """Analyse every topology present. Returns {topology: comp_results}."""
    topos = [t["topology"] for t in config.VALIDATION_TARGETS
             if any(r.get("topology") == t["topology"] for r in records)]
    out = {}
    for topo in topos:
        out[topo] = analyze_topology(records, results_dir, topo)
    return out


# Backward-compatible alias used by older callers.
def analyze_ffc(records: list, results_dir: str) -> dict:
    return analyze_topology(records, results_dir, "ffc")


def write_csv(records: list, results_dir: str):
    """Write full results CSV."""
    csv_path = os.path.join(results_dir, "dft_vs_symmetry.csv")
    fieldnames = ["id", "topology", "composition_label", "composition", "sym_class",
                  "rho", "stabilizer_size", "gamma_size",
                  "dft_e_per_atom", "dft_e_hartree", "n_electrons",
                  "multiplicity", "n_atoms", "probe_cu", "status"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    print(f"[analyze] CSV written to {csv_path}")
