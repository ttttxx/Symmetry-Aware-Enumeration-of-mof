import os, sys, json, random
from collections import Counter, defaultdict
sys.path.insert(0, '.')

# ── NetworkX v3 compat: shortest_path returns generator in v3 ─────────────
import networkx as nx
_nx_sp = nx.shortest_path
nx.shortest_path = lambda *a, **kw: dict(_nx_sp(*a, **kw))

from main import load_graph, enrich_graph_data_with_topology, calculate_topology_requirements
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator
from assembler import StructureAssembler
from dft_validation import config


# ═══════════════════════════════════════════════════════════════════════════
# Pure helpers (unit-tested independently of pymatgen / assembler)
# ═══════════════════════════════════════════════════════════════════════════
def _graph_path(topology: str) -> str:
    return os.path.join("data", topology, f"{topology}_graph.json")


def _apply_cd_substitution(atom_choices: dict) -> dict:
    """Replace 4c_Cu_1_Ch → 4c_Cd_1_Ch in atom_choices if USE_CLOSED_SHELL_VARIANT."""
    if not config.USE_CLOSED_SHELL_VARIANT:
        return atom_choices
    return {
        deg: ["4c_Cd_1_Ch" if l == "4c_Cu_1_Ch" else l for l in ligs]
        for deg, ligs in atom_choices.items()
    }


def comp_key(coloring: list) -> dict:
    """Return {ligand_name: count} composition dict for a coloring."""
    return dict(Counter(coloring))


def comp_label_key(coloring: list) -> str:
    """Canonical, sortable string key for a composition (chemistry fingerprint)."""
    return ",".join(f"{k}:{v}" for k, v in sorted(comp_key(coloring).items()))


def classify_sym_classes(results: list) -> list:
    """
    Assign a GLOBAL sym_class to each result (top-2 stab → High, interior >1 → Med,
    stab==1 → Low). Identical bucketing to the validated ffc path; sym_class is used
    only for plot colouring and the 'lowest-E in High' qualitative check — the
    Spearman x-axis is the continuous rho.
    """
    all_stabs = sorted(set(r["stabilizer_size"] for r in results), reverse=True)
    high_vals, med_vals = set(), set()
    if len(all_stabs) <= 2:
        if all_stabs and all_stabs[0] > 1:
            high_vals.add(all_stabs[0])
    elif len(all_stabs) == 3:
        high_vals.add(all_stabs[0])
        if all_stabs[1] > 1:
            med_vals.add(all_stabs[1])
    else:
        high_vals.update(all_stabs[:2])
        for v in all_stabs[2:]:
            if v > 1:
                med_vals.add(v)

    out = []
    for r in results:
        s = r["stabilizer_size"]
        cls = "High" if s in high_vals else ("Med" if s in med_vals else "Low")
        out.append({**r, "sym_class": cls})
    return out


def group_by_composition(classified: list) -> dict:
    """{composition_key_string: [records]}."""
    groups = defaultdict(list)
    for r in classified:
        groups[comp_label_key(r["coloring"])].append(r)
    return dict(groups)


def _distinct_rho(recs: list) -> int:
    return len(set(r["stabilizer_size"] for r in recs))


def _has_high(recs: list) -> bool:
    return any(r["sym_class"] == "High" for r in recs)


def select_compositions(groups: dict, min_rho_values: int,
                        n_compositions: int, exhaustive: bool) -> list:
    """
    Keep compositions whose INTERNAL rho diversity ≥ min_rho_values (controls the
    chemistry confound: within a composition the Cd:B ratio is fixed, so a rho–E
    trend is a pure symmetry effect). Rank by (rho-richness, contains-High, size).
    Returns ordered list of (comp_key_string, [records]).
    """
    elig = [(k, v) for k, v in groups.items() if _distinct_rho(v) >= min_rho_values]
    elig.sort(key=lambda kv: (-_distinct_rho(kv[1]), not _has_high(kv[1]),
                              -len(kv[1]), kv[0]))
    return elig if exhaustive else elig[:n_compositions]


def stratified_sample_by_rho(recs: list, n_per_rho: int, seed: int) -> list:
    """Sample up to n_per_rho structures from each distinct stabilizer_size (=rho)."""
    rng = random.Random(seed)
    by_stab = defaultdict(list)
    for r in recs:
        by_stab[r["stabilizer_size"]].append(r)
    out = []
    for stab in sorted(by_stab, reverse=True):
        bucket = sorted(by_stab[stab], key=lambda r: r["id"])
        rng.shuffle(bucket)
        out.extend(bucket[:n_per_rho])
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Enumeration + assembly
# ═══════════════════════════════════════════════════════════════════════════
def _enumerate_topology(topology: str, atom_choices: dict):
    """Return (results, graph_data, perms) for a topology + atom_choices."""
    graph_data = enrich_graph_data_with_topology(load_graph(_graph_path(topology)))
    topo_req   = calculate_topology_requirements(graph_data)
    node_degs  = [n.get("degree", 0) for n in graph_data["nodes"]]
    perms      = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
    all_ligs   = sorted(set(l for ls in atom_choices.values() for l in ls))
    results    = MOFEnumerator(
        perms, graph_data["n_nodes"], verbose=False,
        topology_requirements=topo_req,
        user_selected_ligands=atom_choices,
        node_degrees=node_degs,
    ).enumerate(all_ligs)
    return results, graph_data, perms


def _topology_selection_params(topology: str) -> dict:
    """Per-topology knobs, with safe fallbacks if config lacks the new fields."""
    if topology == "ffc":
        return dict(
            exhaustive=getattr(config, "FFC_EXHAUSTIVE", True),
            min_rho=getattr(config, "FFC_MIN_RHO", 2),
            n_comp=10**9,
            n_per_rho=getattr(config, "FFC_N_PER_RHO", 999),
        )
    # fjh (and any other large topology) — sample, don't exhaust
    return dict(
        exhaustive=False,
        min_rho=getattr(config, "FJH_MIN_RHO", 3),
        n_comp=getattr(config, "FJH_N_COMPOSITIONS", 1),
        n_per_rho=getattr(config, "FJH_N_PER_RHO", 2),
    )


def _make_label(idx: int) -> str:
    return f"comp_{chr(ord('a') + idx)}" if idx < 26 else f"comp_{idx}"


def _build_records_for_topology(target: dict):
    """
    Enumerate one topology, pick rho-rich compositions, (sample if not exhaustive),
    and return pending records (pre-assembly) + (graph_data, perms) for reuse.
    """
    topology = target["topology"]
    ac = _apply_cd_substitution(target["atom_choices"])
    results, graph_data, perms = _enumerate_topology(topology, ac)
    gamma = len(perms)
    classified = classify_sym_classes(results)
    groups = group_by_composition(classified)
    p = _topology_selection_params(topology)

    selected = select_compositions(groups, p["min_rho"], p["n_comp"], p["exhaustive"])

    mode_str = "exhaustive" if p["exhaustive"] else \
        f"sample {p['n_per_rho']}/rho x {p['n_comp']} comp"
    print(f"\n[select:{topology}] |Gamma|={gamma}  total_structures={len(results)}  "
          f"compositions={len(groups)}")
    print(f"[select:{topology}] eligible (distinct-rho>={p['min_rho']}): {len(selected)}  "
          f"mode={mode_str}")

    pending = []
    for idx, (ckey, recs) in enumerate(selected):
        label = _make_label(idx)
        chosen = recs if p["exhaustive"] else stratified_sample_by_rho(
            recs, p["n_per_rho"], config.RANDOM_SEED + idx)
        rho_dist = Counter(r["stabilizer_size"] for r in chosen)
        print(f"    {label}: {ckey}")
        print(f"        distinct_rho={_distinct_rho(recs)}  has_High={_has_high(recs)}  "
              f"selected={len(chosen)}  rho(stab) dist={dict(sorted(rho_dist.items(), reverse=True))}")
        for r in chosen:
            pending.append({
                "id": f"{topology}_{label}_{r['id']}",
                "topology": topology,
                "composition_label": label,
                "composition": ckey,
                "coloring": r["coloring"],
                "stabilizer_size": r["stabilizer_size"],
                "gamma_size": gamma,
                "rho": r["stabilizer_size"] / gamma,
                "sym_class": r["sym_class"],
                "n_atoms": None,
                "init_cif": None,
                "probe_cu": False,
                "status": "selected",
            })
    return pending, graph_data, perms


def _assemble_records_with_graph(records: list, graph_data: dict, forced_edge: str,
                                 out_dir: str, topology: str) -> list:
    """Assemble CIFs, update n_atoms + init_cif. Returns updated records."""
    init_dir = os.path.join(out_dir, topology, "init")
    os.makedirs(init_dir, exist_ok=True)
    assembler = StructureAssembler(output_dir=init_dir)

    from pymatgen.core import Structure
    max_atoms = getattr(config, "MAX_ATOMS", 10**9)

    updated = []
    for rec in records:
        smap = {str(i): c for i, c in enumerate(rec["coloring"])}
        cif_name = f"{rec['id']}_rho{rec['rho']:.4f}_{rec['sym_class']}.cif"
        cif_path = os.path.join(init_dir, cif_name)
        try:
            assembler.assemble_and_save(cif_path, smap, graph_data, forced_edge=forced_edge)
            struct = Structure.from_file(cif_path)
            n = len(struct)
            if n > max_atoms:
                print(f"  [select] WARNING {rec['id']}: n_atoms={n} > MAX_ATOMS={max_atoms} "
                      f"(kept; raise MAX_ATOMS or expect a long/heavy DFT job)")
            rec = {**rec, "n_atoms": n, "init_cif": cif_path, "status": "assembled"}
        except Exception as e:
            print(f"  [select] Assembly failed for {rec['id']}: {e}")
            rec = {**rec, "status": f"assembly_failed: {e}"}
        updated.append(rec)
    return updated


def _cu_probe_record(target: dict, ffc_perms) -> dict | None:
    """Pick 1 High-sym ffc structure using ORIGINAL Cu atom_choices (open-shell)."""
    if not config.CU_PROBE:
        return None
    ac = target["atom_choices"]  # no Cd substitution
    gamma = len(ffc_perms)
    results, _, _ = _enumerate_topology(target["topology"], ac)
    classified = classify_sym_classes(results)
    high = [r for r in classified if r["sym_class"] == "High"]
    if not high:
        return None
    best = max(high, key=lambda r: (r["stabilizer_size"], -r["id"]))
    return {
        "id": "ffc_cu_probe_0",
        "topology": "ffc",
        "composition_label": "cu_probe",
        "composition": comp_label_key(best["coloring"]),
        "coloring": best["coloring"],
        "stabilizer_size": best["stabilizer_size"],
        "gamma_size": gamma,
        "rho": best["stabilizer_size"] / gamma,
        "sym_class": best["sym_class"],
        "n_atoms": None,
        "init_cif": None,
        "probe_cu": True,
        "status": "selected",
    }


# ═══════════════════════════════════════════════════════════════════════════
# Main entry
# ═══════════════════════════════════════════════════════════════════════════
def run_selection(results_dir: str = None) -> list:
    """
    Phase 1. Iterate VALIDATION_TARGETS in order (ffc first, fast → fjh last).
    For each topology select rho-rich compositions, assemble, accumulate into ONE
    selected_structures.json. ffc Cu probe appended last.
    """
    if results_dir is None:
        results_dir = config.RESULTS_DIR
    os.makedirs(results_dir, exist_ok=True)

    all_records = []
    ffc_perms = None

    for target in config.VALIDATION_TARGETS:
        topology = target["topology"]
        pending, graph_data, perms = _build_records_for_topology(target)
        if topology == "ffc":
            ffc_perms = perms
        assembled = _assemble_records_with_graph(
            pending, graph_data, target["forced_edge"], results_dir, topology)
        all_records.extend(assembled)

    # ffc Cu probe (open-shell) — re-enumerate with original Cu choices
    ffc_target = next((t for t in config.VALIDATION_TARGETS if t["topology"] == "ffc"), None)
    if ffc_target is not None and ffc_perms is not None:
        probe = _cu_probe_record(ffc_target, ffc_perms)
        if probe:
            graph_data = enrich_graph_data_with_topology(load_graph(_graph_path("ffc")))
            probe_assembled = _assemble_records_with_graph(
                [probe], graph_data, ffc_target["forced_edge"], results_dir, "ffc")
            all_records.extend(probe_assembled)

    json_path = os.path.join(results_dir, "selected_structures.json")
    with open(json_path, "w") as f:
        json.dump(all_records, f, indent=2)

    assembled = sum(1 for r in all_records if r["status"] == "assembled")
    print(f"\n[select] Wrote {len(all_records)} records to {json_path}")
    print(f"[select] Assembled OK: {assembled}/{len(all_records)}")
    by_topo = Counter(r["topology"] for r in all_records)
    print(f"[select] By topology: {dict(by_topo)}")
    return all_records
