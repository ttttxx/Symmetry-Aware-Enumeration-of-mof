import os, sys, json, random
from collections import Counter
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


def _comp_matches(coloring: list, target: dict) -> bool:
    return comp_key(coloring) == target


def classify_sym_classes(results: list) -> list:
    """
    Assign sym_class to each result using physical_validation.py bucketing logic
    (lines 105-120): top-2 stab values → High, interior >1 → Med, stab==1 → Low.
    Returns new list with 'sym_class' key added.
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
        if s in high_vals:
            cls = "High"
        elif s in med_vals:
            cls = "Med"
        else:
            cls = "Low"
        out.append({**r, "sym_class": cls})
    return out


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


def enumerate_ffc(target: dict) -> list:
    """
    Enumerate ffc, apply Cd substitution, filter to FFC_COMPOSITIONS,
    classify sym_classes, return all matching records (no sampling).
    """
    ac = _apply_cd_substitution(target["atom_choices"])
    results, graph_data, perms = _enumerate_topology(target["topology"], ac)
    gamma_size = len(perms)
    classified = classify_sym_classes(results)

    records = []
    for comp_dict, label in zip(config.FFC_COMPOSITIONS, config.FFC_COMP_LABELS):
        matching = [r for r in classified if _comp_matches(r["coloring"], comp_dict)]
        for r in matching:
            records.append({
                "id": f"ffc_{label}_{r['id']}",
                "topology": "ffc",
                "composition_label": label,
                "coloring": r["coloring"],
                "stabilizer_size": r["stabilizer_size"],
                "gamma_size": gamma_size,
                "rho": r["stabilizer_size"] / gamma_size,
                "sym_class": r["sym_class"],
                "n_atoms": None,  # filled after assembly
                "init_cif": None,
                "probe_cu": False,
                "status": "selected",
            })
    return records


def _cu_probe_record(target: dict, ffc_graph_data, ffc_perms) -> dict | None:
    """Pick 1 High-sym ffc structure using original Cu atom_choices."""
    if not config.CU_PROBE:
        return None
    # Use ORIGINAL atom_choices (no Cd substitution)
    ac = target["atom_choices"]
    gamma_size = len(ffc_perms)
    results, _, _ = _enumerate_topology(target["topology"], ac)
    classified = classify_sym_classes(results)
    high = [r for r in classified if r["sym_class"] == "High"]
    if not high:
        return None
    # Pick highest stabilizer_size first, then by id for determinism
    best = max(high, key=lambda r: (r["stabilizer_size"], -r["id"]))
    return {
        "id": "ffc_cu_probe_0",
        "topology": "ffc",
        "composition_label": "cu_probe",
        "coloring": best["coloring"],
        "stabilizer_size": best["stabilizer_size"],
        "gamma_size": gamma_size,
        "rho": best["stabilizer_size"] / gamma_size,
        "sym_class": best["sym_class"],
        "n_atoms": None,
        "init_cif": None,
        "probe_cu": True,
        "status": "selected",
    }


def _assemble_records_with_graph(records: list, graph_data: dict, forced_edge: str,
                                  out_dir: str, topology: str) -> list:
    """Assemble CIFs, update n_atoms + init_cif. Returns updated records."""
    init_dir = os.path.join(out_dir, topology, "init")
    os.makedirs(init_dir, exist_ok=True)
    assembler = StructureAssembler(output_dir=init_dir)

    from pymatgen.core import Structure

    updated = []
    for rec in records:
        smap = {str(i): c for i, c in enumerate(rec["coloring"])}
        cif_name = f"{rec['id']}_rho{rec['rho']:.4f}_{rec['sym_class']}.cif"
        cif_path = os.path.join(init_dir, cif_name)
        try:
            assembler.assemble_and_save(cif_path, smap, graph_data, forced_edge=forced_edge)
            struct = Structure.from_file(cif_path)
            rec = {**rec, "n_atoms": len(struct), "init_cif": cif_path, "status": "assembled"}
        except Exception as e:
            print(f"  [select] Assembly failed for {rec['id']}: {e}")
            rec = {**rec, "status": f"assembly_failed: {e}"}
        updated.append(rec)
    return updated


def run_selection(results_dir: str = None) -> list:
    """
    Main entry point for Phase 1.
    Enumerates ffc Comp A+B + Cu probe, assembles CIFs, writes selected_structures.json.
    Returns the list of records.
    """
    if results_dir is None:
        results_dir = config.RESULTS_DIR
    os.makedirs(results_dir, exist_ok=True)

    all_records = []

    ffc_target = next(t for t in config.VALIDATION_TARGETS if t["topology"] == "ffc")
    ac_cd = _apply_cd_substitution(ffc_target["atom_choices"])

    # Enumerate once (Cd variant) — reuse graph_data and perms for probe too
    graph_data = enrich_graph_data_with_topology(load_graph(_graph_path("ffc")))
    topo_req   = calculate_topology_requirements(graph_data)
    node_degs  = [n.get("degree", 0) for n in graph_data["nodes"]]
    perms      = SymmetryEngine(graph_data, verbose=False).compute_automorphism_group()
    gamma_size = len(perms)
    all_ligs   = sorted(set(l for ls in ac_cd.values() for l in ls))
    results    = MOFEnumerator(
        perms, graph_data["n_nodes"], verbose=False,
        topology_requirements=topo_req,
        user_selected_ligands=ac_cd,
        node_degrees=node_degs,
    ).enumerate(all_ligs)
    classified = classify_sym_classes(results)

    cd_pending = []
    for comp_dict, label in zip(config.FFC_COMPOSITIONS, config.FFC_COMP_LABELS):
        matching = [r for r in classified if _comp_matches(r["coloring"], comp_dict)]
        for r in matching:
            cd_pending.append({
                "id": f"ffc_{label}_{r['id']}",
                "topology": "ffc",
                "composition_label": label,
                "coloring": r["coloring"],
                "stabilizer_size": r["stabilizer_size"],
                "gamma_size": gamma_size,
                "rho": r["stabilizer_size"] / gamma_size,
                "sym_class": r["sym_class"],
                "n_atoms": None,
                "init_cif": None,
                "probe_cu": False,
                "status": "selected",
            })

    # Assemble all Cd records
    cd_records = _assemble_records_with_graph(
        cd_pending, graph_data, ffc_target["forced_edge"], results_dir, "ffc"
    )
    all_records.extend(cd_records)

    # Cu probe (re-enumerate with original Cu choices)
    probe = _cu_probe_record(ffc_target, graph_data, perms)
    if probe:
        probe_records = _assemble_records_with_graph(
            [probe], graph_data, ffc_target["forced_edge"], results_dir, "ffc"
        )
        all_records.extend(probe_records)

    json_path = os.path.join(results_dir, "selected_structures.json")
    with open(json_path, "w") as f:
        json.dump(all_records, f, indent=2)
    print(f"[select] Wrote {len(all_records)} records to {json_path}")
    assembled = sum(1 for r in all_records if r["status"] == "assembled")
    print(f"[select] Assembled: {assembled}/{len(all_records)}")
    return all_records
