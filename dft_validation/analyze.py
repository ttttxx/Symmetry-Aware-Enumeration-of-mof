import os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from dft_validation import config


SYM_COLORS = {"High": "#d62728", "Med": "#1f77b4", "Low": "#95a5a6"}


def _spearman_for_comp(records: list) -> tuple:
    """
    Compute Spearman rs for (rho, dft_e_per_atom).
    Returns (rs, p_val, lowest_in_high: bool).
    """
    converged = [r for r in records if r.get("dft_e_per_atom") is not None]
    if len(converged) < 2:
        return (float("nan"), float("nan"), False)
    rhos     = [r["rho"] for r in converged]
    energies = [r["dft_e_per_atom"] for r in converged]
    rs, p    = spearmanr(rhos, energies)
    min_rec  = min(converged, key=lambda r: r["dft_e_per_atom"])
    return (rs, p, min_rec["sym_class"] == "High")


def get_fjh_gate(comp_results: dict) -> tuple:
    """
    Qualitative gate: both comps must have lowest-E in High AND rs < 0.
    Returns (passed: bool, reason: str).
    """
    for label, stats in comp_results.items():
        if not stats.get("lowest_in_high", False):
            return False, f"{label}: lowest-energy structure is NOT in High class"
        if stats.get("spearman_rs", 0) >= 0:
            return False, (f"{label}: Spearman rs={stats['spearman_rs']:.3f} >= 0 "
                           f"(energy does not decrease with symmetry)")
    return True, "both compositions pass gate"


def analyze_ffc(records: list, results_dir: str) -> dict:
    """
    Per-composition Spearman + side-by-side rho-E plot.
    Returns dict comp_label → {spearman_rs, p_val, lowest_in_high}.
    """
    comp_results = {}
    labels = sorted({r["composition_label"] for r in records
                     if r.get("topology") == "ffc" and not r.get("probe_cu")})

    fig, axes = plt.subplots(1, len(labels), figsize=(6 * len(labels), 5), sharey=False)
    if len(labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, labels):
        subset = [r for r in records
                  if r.get("composition_label") == label
                  and r.get("dft_e_per_atom") is not None]
        rs, p, lowest_in_high = _spearman_for_comp(subset)
        comp_results[label] = {"spearman_rs": rs, "p_val": p, "lowest_in_high": lowest_in_high}

        for cls in ["High", "Med", "Low"]:
            pts = [r for r in subset if r.get("sym_class") == cls]
            if pts:
                ax.scatter([r["rho"] for r in pts],
                           [r["dft_e_per_atom"] for r in pts],
                           color=SYM_COLORS[cls], label=cls, s=60, alpha=0.85, zorder=3)

        ax.set_xlabel(r"Symmetry Index $\rho$", fontsize=12)
        ax.set_ylabel("DFT Energy (eV/atom)", fontsize=12)
        rs_str = f"{rs:.3f}" if rs == rs else "n/a"  # NaN check
        ax.set_title(f"ffc {label}  (Spearman $\\rho_s$={rs_str}, p={p:.3f})", fontsize=11)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(results_dir, "rho_energy_ffc.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[analyze] Figure saved to {fig_path}")
    return comp_results


def write_csv(records: list, results_dir: str):
    """Write full results CSV."""
    csv_path = os.path.join(results_dir, "dft_vs_symmetry.csv")
    fieldnames = ["id", "topology", "composition_label", "sym_class", "rho",
                  "dft_e_per_atom", "dft_e_hartree", "n_electrons",
                  "multiplicity", "n_atoms", "status"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    print(f"[analyze] CSV written to {csv_path}")
