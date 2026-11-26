# graph/build_graph.py
# -*- coding: utf-8 -*-
"""
Step2: Build attribute graph from Step1 output (structure json)
Saves: <prefix>_graph.json  (and <prefix>_graph.gpickle)
Optional: saves <prefix>_graph.png if --plot

Usage:
  python utils/build_graph.py data/pcu --plot
"""

import os
import sys
import json
import argparse
import pickle
from typing import Tuple, Dict, Any, List, Set
import numpy as np

# pymatgen / networkx / matplotlib
from pymatgen.core import Structure, Lattice
from pymatgen.analysis.local_env import CrystalNN
import networkx as nx
import matplotlib.pyplot as plt

# ---------- Helpers ----------
ORGANIC_SET = {"C", "H", "O", "N", "S", "P", "B", "Cl", "F", "Br", "I"}


def find_prefix_and_files(path: str) -> Tuple[str, str]:
    path = os.path.abspath(path)
    if os.path.isdir(path):
        for fname in os.listdir(path):
            if "structure" in fname and fname.lower().endswith(".json"):
                return path, os.path.join(path, fname)
        for fname in os.listdir(path):
            if fname.lower().endswith("_topology.json"):
                prefix = fname[:-len("_topology.json")]
                sfile = os.path.join(path, prefix + "_structure_info.json")
                if os.path.exists(sfile):
                    return path, sfile
        raise FileNotFoundError(f"No structure json found in directory: {path}")
    else:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input path does not exist: {path}")
        dirname = os.path.dirname(path)
        basename = os.path.basename(path)
        if "structure" in basename.lower() and basename.lower().endswith(".json"):
            return dirname, path
        if basename.lower().endswith("_topology.json"):
            prefix = basename[:-len("_topology.json")]
            sfile = os.path.join(dirname, prefix + "_structure_info.json")
            if os.path.exists(sfile):
                return dirname, sfile
        if basename.lower().endswith(".json"):
            return dirname, path
        raise FileNotFoundError("Unsupported input file.")


def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def rebuild_structure_from_structure_json(struct_json: Dict[str, Any]) -> Structure:
    lat = struct_json.get("lattice_matrix")
    if lat is None:
        raise ValueError("structure json missing 'lattice_matrix'")
    lattice = Lattice(lat)
    species = struct_json.get("species")
    frac_coords = struct_json.get("frac_coords")
    if (species is None or frac_coords is None) and "nodes" in struct_json:
        nodes = struct_json["nodes"]
        species = [n["element"] for n in nodes]
        frac_coords = [n.get("frac") or n.get("coord") for n in nodes]
    if species is None or frac_coords is None:
        raise ValueError("structure json must contain 'species' and 'frac_coords'")
    return Structure(lattice, species, frac_coords)


def node_type_from_element(el: str) -> str:
    if el is None: return "unknown"
    # Special case for dummy atoms often found in topologies
    if "Z" in el and "+" in el: return "metal"  # Assume topology nodes are metal nodes usually

    e = el.capitalize()
    return "linker" if e in ORGANIC_SET else "metal"


# ---------- Connection Logic ----------

def _get_edges_crystalnn(structure: Structure) -> List[Dict]:
    """Use Pymatgen CrystalNN to find chemical bonds."""
    cnn = CrystalNN()
    edges = []
    added = set()
    N = len(structure.sites)

    for i in range(N):
        try:
            # CrystalNN might fail on dummy elements
            neighbors = cnn.get_nn_info(structure, i)
        except Exception:
            continue

        for nn in neighbors:
            j = nn.get("site_index")
            if j is None or i == j: continue

            # Sort to handle undirected edges unique
            u, v = sorted((int(i), int(j)))
            if (u, v) in added: continue

            try:
                dist = structure.get_distance(i, j)
            except:
                dist = 0.0

            edges.append({
                "source": u,
                "target": v,
                "length": float(dist),
                "bond_type": "coordination"
            })
            added.add((u, v))
    return edges


def _get_edges_by_distance(structure: Structure, tolerance: float = 1.10) -> List[Dict]:
    """
    Fallback method: Connect nodes to their 'first shell' neighbors.
    Finds the nearest neighbor distance (d_min) for each node,
    and connects all neighbors within d_min * tolerance.
    """
    edges = []
    added = set()
    N = len(structure.sites)

    print(f"  [Info] Using Distance-based connection (tolerance={tolerance})...")

    for i in range(N):
        # Get neighbors within a generous radius (e.g. 20 Angstrom)
        # We assume topology edges aren't longer than 20A
        neighbors = structure.get_neighbors(structure[i], r=20.0)

        if not neighbors:
            continue

        # Sort by distance
        neighbors.sort(key=lambda x: x.nn_distance)

        # Determine local cutoff based on the closest neighbor
        # Use first neighbor > 0.1A (to avoid self-image issues if any)
        closest_dist = neighbors[0].nn_distance
        cutoff = closest_dist * tolerance

        for n in neighbors:
            d = n.nn_distance
            if d > cutoff:
                break  # Moved past the first shell

            j = n.index
            if i == j: continue

            u, v = sorted((int(i), int(j)))
            if (u, v) in added: continue

            edges.append({
                "source": u,
                "target": v,
                "length": float(d),
                "bond_type": "topology_link"
            })
            added.add((u, v))

    return edges


# ---------- Main Graph Builder ----------
def build_graph_from_structure(structure: Structure,
                               only_metal_linker: bool = True,
                               method: str = "crystalnn") -> Dict[str, Any]:
    N = len(structure.sites)
    nodes = []
    for i, site in enumerate(structure.sites):
        el = str(site.specie)
        frac = [float(x) for x in site.frac_coords]
        ntype = node_type_from_element(el)
        nodes.append({"id": int(i), "element": el, "coord": frac, "type": ntype, "degree": 0})

    # 1. Try specified method (default CrystalNN)
    edges = []
    if method == "crystalnn":
        edges = _get_edges_crystalnn(structure)

    # 2. Fallback: If CrystalNN returned empty (common for Z0+ dummy atoms), use Distance
    if len(edges) == 0:
        if method == "crystalnn":
            print("  [Warning] CrystalNN found 0 edges (dummy elements?). Switching to geometry-based connection.")
        edges = _get_edges_by_distance(structure)

    # 3. Optional Filter: Metal-Linker Only
    # For pure topologies (all Z0+), we should SKIP this filter, otherwise everything is deleted.
    # We check if we have mixed types.
    types_present = set(n["type"] for n in nodes)
    has_metal_and_linker = ("metal" in types_present and "linker" in types_present)

    final_edges = []
    if only_metal_linker and has_metal_and_linker:
        # standard MOF case
        for e in edges:
            t1 = nodes[e["source"]]["type"]
            t2 = nodes[e["target"]]["type"]
            # keep only metal-linker
            if {t1, t2} == {"metal", "linker"}:
                final_edges.append(e)
    else:
        # For pure topologies (all "metal" or all "unknown"), keep all edges
        final_edges = edges

    # Update degrees
    deg = {i: 0 for i in range(N)}
    for e in final_edges:
        deg[e["source"]] += 1
        deg[e["target"]] += 1
    for node in nodes:
        node["degree"] = deg.get(node["id"], 0)

    graph = {"n_nodes": N, "nodes": nodes, "edges": final_edges}
    return graph


def save_graph_json(graph: Dict[str, Any], out_json_path: str):
    with open(out_json_path, "w", encoding="utf-8") as f:
        json.dump(graph, f, indent=2, ensure_ascii=False)


def save_graph_gpickle(graph: Dict[str, Any], out_gpickle_path: str):
    G = nx.Graph()
    for n in graph["nodes"]:
        G.add_node(n["id"], **{k: v for k, v in n.items() if k != "id"})
    for e in graph["edges"]:
        G.add_edge(e["source"], e["target"], **{k: v for k, v in e.items() if k not in ("source", "target")})
    with open(out_gpickle_path, 'wb') as f:
        pickle.dump(G, f)


def plot_graph(graph: Dict[str, Any], out_png: str, figsize=(8, 8), show_labels=True):
    G = nx.Graph()
    for n in graph["nodes"]:
        G.add_node(n["id"], element=n["element"], type=n["type"])
    for e in graph["edges"]:
        G.add_edge(e["source"], e["target"])

    pos = nx.spring_layout(G, seed=42)
    types = nx.get_node_attributes(G, "type")
    color_map = []
    for node in G.nodes():
        t = types.get(node, "unknown")
        if t == "metal":
            color_map.append("#ff7f0e")
        elif t == "linker":
            color_map.append("#1f77b4")
        else:
            color_map.append("#7f7f7f")

    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(G, pos, node_color=color_map, node_size=300)
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.8)
    if show_labels:
        labels = {n: str(n) for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=8)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build attribute graph from structure json")
    parser.add_argument("input", help="prefix dir or structure json")
    parser.add_argument("--out_root", default=None)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--method", choices=["crystalnn", "distance"], default="crystalnn")
    parser.add_argument("--no_metal_linker", action="store_true")
    args = parser.parse_args(argv)

    try:
        prefix_dir, struct_json_path = find_prefix_and_files(args.input)
    except Exception as e:
        print("Error locating structure json:", e)
        sys.exit(2)

    out_dir = args.out_root if args.out_root else prefix_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    prefix_name = os.path.basename(os.path.normpath(prefix_dir))

    struct_json = load_json(struct_json_path)
    structure = rebuild_structure_from_structure_json(struct_json)

    print(f"Building graph for '{prefix_name}'...")

    # Auto-detect if we should force no_metal_linker (e.g. if it's a raw topology)
    # The function build_graph_from_structure now handles the fallback logic internally
    graph = build_graph_from_structure(structure,
                                       only_metal_linker=not args.no_metal_linker,
                                       method=args.method)

    print(f"  Nodes: {graph['n_nodes']}, Edges: {len(graph['edges'])}")

    graph_json_path = os.path.join(out_dir, f"{prefix_name}_graph.json")
    graph_gpickle_path = os.path.join(out_dir, f"{prefix_name}_graph.gpickle")

    save_graph_json(graph, graph_json_path)
    save_graph_gpickle(graph, graph_gpickle_path)
    print(f"Wrote: {graph_json_path}")

    if args.plot:
        png_path = os.path.join(out_dir, f"{prefix_name}_graph.png")
        try:
            plot_graph(graph, png_path)
            print(f"Wrote plot: {png_path}")
        except Exception as e:
            print("Plot failed:", e)


if __name__ == "__main__":
    main()