import os, csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
from dft_validation import config

RHO_CMAP = "viridis"
RHO_BANDS = [
    ("Low", 0.0, 0.25, "#D8D8D8"),
    ("Medium", 0.25, 0.75, "#DCEBF7"),
    ("High", 0.75, 1.0, "#E3F1DE"),
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

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
    Qualitative retention gate (kept for backward compat / informational use):
    composition-wise minima should be retained outside the low-rho interval.
    """
    if not comp_results:
        return None, "no compositions analysed"
    for label, stats in comp_results.items():
        if stats.get("n", 0) < 2:
            return None, f"{label}: insufficient converged DFT points"
        if stats.get("lowest_band") == "Low":
            rho = stats.get("lowest_rho")
            return False, f"{label}: lowest-energy structure falls in Low rho (rho={rho:g})"
    return True, "all composition-wise minima fall in Medium/High rho"

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
        min_rec = min(subset, key=lambda r: r["dft_e_per_atom"]) if subset else None
        comp_results[label] = {"spearman_rs": rs, "p_val": p,
                               "lowest_at_max_rho": lowest_at_max_rho,
                               "lowest_rho": min_rec.get("rho") if min_rec else None,
                               "lowest_band": _rho_band(min_rec["rho"]) if min_rec else None,
                               "n": n}

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

    if topology in {"ffc", "fjh"}:
        _plot_pooled_topology(records, results_dir, topology)
    return comp_results

def _rho_band(rho: float) -> str:
    if rho < 0.25:
        return "Low"
    if rho <= 0.75:
        return "Medium"
    return "High"

def _add_rho_bands(ax) -> None:
    for label, lo, hi, color in RHO_BANDS:
        ax.axvspan(lo, hi, color=color, alpha=0.42, lw=0, zorder=0)
        ax.text((lo + hi) / 2, 0.985, label, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.3, color="#4D4D4D")

def _composition_minima(records: list) -> list:
    groups = defaultdict(list)
    for rec in records:
        groups[rec.get("composition_label", "")].append(rec)
    return [min(group, key=lambda r: r["dft_e_per_atom"])
            for group in groups.values() if group]

def _normalized_records(records: list) -> list:
    groups = defaultdict(list)
    for rec in records:
        groups[rec.get("composition_label", "")].append(rec)
    normalized = []
    for group in groups.values():
        baseline = min(r["dft_e_per_atom"] for r in group)
        for rec in group:
            normalized.append({**rec, "delta_mev_atom":
                               (rec["dft_e_per_atom"] - baseline) * 1000.0})
    return normalized

def _plot_single_screening_panel(ax, points: list, minima: list, *,
                                 y_key: str, y_label: str,
                                 legend: bool = False) -> None:
    colors = {"all": "#5F6F7A", "min": "#B64342"}
    comp_labels = sorted({r.get("composition_label", "") for r in points})
    marker_map = {label: marker for label, marker in zip(
        comp_labels, ["o", "s", "^", "D", "P", "X", "v", "<", ">", "h", "*"]
    )}

    _add_rho_bands(ax)
    ax.set_xlim(0.04, 1.06)
    ax.set_xticks([0.125, 0.25, 0.5, 1.0])
    ax.set_xlabel(r"Symmetry index, $\rho$")
    ax.set_ylabel(y_label)
    ax.grid(axis="y", color="#D8D8D8", linewidth=0.5, alpha=0.75, zorder=0)

    for label in comp_labels:
        subset = [r for r in points if r.get("composition_label") == label]
        ax.scatter([r["rho"] for r in subset],
                   [r[y_key] for r in subset],
                   marker=marker_map[label], s=30, color=colors["all"],
                   alpha=0.58, edgecolor="white", linewidth=0.35, zorder=2)

    ax.scatter([r["rho"] for r in minima],
               [r[y_key] for r in minima],
               s=68, facecolor="white", edgecolor=colors["min"],
               linewidth=1.3, zorder=4)

    if legend:
        legend_handles = [
            Line2D([0], [0], marker="o", linestyle="none", markersize=4.8,
                   markerfacecolor=colors["all"], markeredgecolor="white",
                   alpha=0.65, label="completed candidate"),
            Line2D([0], [0], marker="o", linestyle="none", markersize=6.8,
                   markerfacecolor="white", markeredgecolor=colors["min"],
                   markeredgewidth=1.25, label="composition minimum"),
        ]
        ax.legend(handles=legend_handles, loc="upper left",
                  bbox_to_anchor=(0.0, -0.2), ncol=2,
                  handletextpad=0.4, columnspacing=1.3)

def _save_publication_figure(fig, path_stem: str) -> None:
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{path_stem}.{ext}", dpi=600 if ext == "png" else None,
                    bbox_inches="tight")

def _plot_single_composition_publication_figure(records: list, results_dir: str,
                                                topology: str, composition_label: str) -> None:
    """Publication-oriented figure for one fully sampled composition."""
    pts = [r for r in records
           if r.get("topology") == topology
           and r.get("composition_label") == composition_label
           and not r.get("probe_cu")
           and r.get("dft_e_per_atom") is not None]
    if not pts:
        return

    minima = _composition_minima(pts)
    normalized = _normalized_records(pts)
    rs, p = spearmanr([r["rho"] for r in pts], [r["dft_e_per_atom"] for r in pts])
    min_rec = min(pts, key=lambda r: r["dft_e_per_atom"])
    rho1 = [r for r in pts if r["rho"] == max(r["rho"] for r in pts)]
    delta_max_rho = ((rho1[0]["dft_e_per_atom"] - min_rec["dft_e_per_atom"]) * 1000.0
                     if rho1 else float("nan"))

    fig, (ax_abs, ax_delta) = plt.subplots(
        1, 2, figsize=(7.2, 3.05), gridspec_kw={"width_ratios": [1.02, 1.0]}
    )
    _plot_single_screening_panel(
        ax_abs, pts, minima,
        y_key="dft_e_per_atom",
        y_label="Single-point DFT energy (eV atom$^{-1}$)",
        legend=True,
    )
    ax_abs.text(0.02, 0.03,
                f"{composition_label}, n={len(pts)}\n"
                f"$\\rho_s$={rs:.2f}, p={p:.2f}",
                transform=ax_abs.transAxes, ha="left", va="bottom",
                fontsize=6.5, color="#4D4D4D")

    minima_norm = [rec for rec in normalized if rec["delta_mev_atom"] == 0]
    _plot_single_screening_panel(
        ax_delta, normalized, minima_norm,
        y_key="delta_mev_atom",
        y_label=r"$\Delta E$ from composition minimum (meV atom$^{-1}$)",
    )
    ax_delta.axhline(0, color="#B64342", linewidth=0.75, alpha=0.65, zorder=1)
    ax_delta.text(0.02, 0.07,
                  "Single-composition minima:\n"
                  f"$\\rho_{{min}}$={min_rec['rho']:g}, "
                  f"$\\Delta E_{{\\rho=1}}$={delta_max_rho:.1f} meV atom$^{{-1}}$",
                  transform=ax_delta.transAxes, ha="left", va="bottom",
                  fontsize=6.5, color="#4D4D4D")
    fig.tight_layout(pad=0.55, w_pad=1.4)

    path_stem = os.path.join(results_dir, f"rho_energy_{topology}_{composition_label}")
    _save_publication_figure(fig, path_stem)
    plt.close(fig)
    print(f"[analyze:{topology}] Single-composition figure saved to {path_stem}.*")

def _plot_fjh_dual_delta_publication_figure(records: list, results_dir: str) -> None:
    """Publication-oriented dual panel delta-E figure for the fully sampled fjh sets."""
    comp_labels = []
    for label in sorted({r.get("composition_label", "") for r in records if r.get("topology") == "fjh"}):
        subset = [r for r in records
                  if r.get("topology") == "fjh"
                  and r.get("composition_label") == label
                  and not r.get("probe_cu")
                  and r.get("dft_e_per_atom") is not None]
        if len(subset) >= 5 and len({r["rho"] for r in subset}) >= 5:
            comp_labels.append(label)

    if len(comp_labels) < 2:
        return

    fig, axes = plt.subplots(1, len(comp_labels), figsize=(7.2, 3.05), sharey=True)
    if len(comp_labels) == 1:
        axes = [axes]

    for idx, (ax, label) in enumerate(zip(axes, comp_labels)):
        pts = [r for r in records
               if r.get("topology") == "fjh"
               and r.get("composition_label") == label
               and not r.get("probe_cu")
               and r.get("dft_e_per_atom") is not None]
        normalized = _normalized_records(pts)
        minima_norm = [rec for rec in normalized if rec["delta_mev_atom"] == 0]
        min_rec = min(pts, key=lambda r: r["dft_e_per_atom"])

        _plot_single_screening_panel(
            ax, normalized, minima_norm,
            y_key="delta_mev_atom",
            y_label=r"$\Delta E$ from composition minimum (meV atom$^{-1}$)",
            legend=(idx == 0),
        )
        ax.axhline(0, color="#B64342", linewidth=0.75, alpha=0.65, zorder=1)
        ax.text(0.02, 0.07,
                f"{label}, n={len(pts)}\n"
                f"$\\rho_{{min}}$={min_rec['rho']:g}",
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize=6.5, color="#4D4D4D")

    fig.tight_layout(pad=0.55, w_pad=1.4)
    path_stem = os.path.join(results_dir, "rho_energy_fjh_dual_delta")
    _save_publication_figure(fig, path_stem)
    plt.close(fig)
    print(f"[analyze:fjh] Dual-composition delta figure saved to {path_stem}.*")

def _plot_pooled_topology(records: list, results_dir: str, topology: str) -> None:
    """Publication-oriented pooled preview with separate absolute and delta panels."""
    pts = [r for r in records
           if r.get("topology") == topology
           and not r.get("probe_cu")
           and r.get("dft_e_per_atom") is not None]
    if not pts:
        return
    minima = _composition_minima(pts)
    normalized = _normalized_records(pts)
    rs, p = spearmanr([r["rho"] for r in pts], [r["dft_e_per_atom"] for r in pts])

    min_band_counts = defaultdict(int)
    for rec in minima:
        min_band_counts[_rho_band(rec["rho"])] += 1
    minima_norm = []
    min_norm_by_label = {
        rec.get("composition_label", ""): rec for rec in normalized if rec["delta_mev_atom"] == 0
    }
    for rec in minima:
        minima_norm.append(min_norm_by_label[rec.get("composition_label", "")])

    fig_abs, ax_abs = plt.subplots(figsize=(3.55, 3.0))
    _plot_single_screening_panel(
        ax_abs, pts, minima,
        y_key="dft_e_per_atom",
        y_label="Single-point DFT energy (eV atom$^{-1}$)",
        legend=True,
    )
    ax_abs.text(0.02, 0.03,
                f"n={len(pts)} candidates, {len(minima)} compositions\n"
                f"pooled $\\rho_s$={rs:.2f}, p={p:.2f}",
                transform=ax_abs.transAxes, ha="left", va="bottom",
                fontsize=6.5, color="#4D4D4D")
    fig_abs.tight_layout(pad=0.45)

    fig_delta, ax_delta = plt.subplots(figsize=(3.55, 3.0))
    _plot_single_screening_panel(
        ax_delta, normalized, minima_norm,
        y_key="delta_mev_atom",
        y_label=r"$\Delta E$ from composition minimum (meV atom$^{-1}$)",
    )
    ax_delta.axhline(0, color="#B64342", linewidth=0.75, alpha=0.65, zorder=1)
    ax_delta.text(0.02, 0.07,
                  "Composition-wise minima:\n"
                  f"Low {min_band_counts['Low']} | "
                  f"Medium {min_band_counts['Medium']} | "
                  f"High {min_band_counts['High']}",
                  transform=ax_delta.transAxes, ha="left", va="bottom",
                  fontsize=6.5, color="#4D4D4D")
    fig_delta.tight_layout(pad=0.45)

    pooled_stem = os.path.join(results_dir, f"rho_energy_{topology}_pooled")
    _save_publication_figure(fig_abs, pooled_stem + "_absolute")
    _save_publication_figure(fig_delta, pooled_stem + "_delta")

    # Keep the legacy combined export as a convenience preview.
    fig_combo, (combo_abs, combo_delta) = plt.subplots(
        1, 2, figsize=(7.2, 3.05), gridspec_kw={"width_ratios": [1.05, 1.0]}
    )
    _plot_single_screening_panel(
        combo_abs, pts, minima,
        y_key="dft_e_per_atom",
        y_label="Single-point DFT energy (eV atom$^{-1}$)",
        legend=True,
    )
    combo_abs.text(0.02, 0.03,
                   f"n={len(pts)} candidates, {len(minima)} compositions\n"
                   f"pooled $\\rho_s$={rs:.2f}, p={p:.2f}",
                   transform=combo_abs.transAxes, ha="left", va="bottom",
                   fontsize=6.5, color="#4D4D4D")
    _plot_single_screening_panel(
        combo_delta, normalized, minima_norm,
        y_key="delta_mev_atom",
        y_label=r"$\Delta E$ from composition minimum (meV atom$^{-1}$)",
    )
    combo_delta.axhline(0, color="#B64342", linewidth=0.75, alpha=0.65, zorder=1)
    combo_delta.text(0.02, 0.07,
                     "Composition-wise minima:\n"
                     f"Low {min_band_counts['Low']} | "
                     f"Medium {min_band_counts['Medium']} | "
                     f"High {min_band_counts['High']}",
                     transform=combo_delta.transAxes, ha="left", va="bottom",
                     fontsize=6.5, color="#4D4D4D")
    fig_combo.tight_layout(pad=0.55, w_pad=1.4)
    _save_publication_figure(fig_combo, pooled_stem)

    plt.close(fig_abs)
    plt.close(fig_delta)
    plt.close(fig_combo)
    print(f"[analyze:{topology}] Pooled descriptive figures saved to {pooled_stem}_*.png")

    if topology == "fjh":
        comp_labels = sorted({r.get("composition_label", "") for r in pts})
        for label in comp_labels:
            subset = [r for r in pts if r.get("composition_label") == label]
            if len(subset) >= 5 and len({r["rho"] for r in subset}) >= 5:
                _plot_single_composition_publication_figure(
                    records, results_dir, topology, label
                )
        _plot_fjh_dual_delta_publication_figure(records, results_dir)

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
