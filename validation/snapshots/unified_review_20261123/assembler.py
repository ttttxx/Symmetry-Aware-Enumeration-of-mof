import os
import sys
import json
import numpy as np
import re
import networkx as nx
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from external_lib.utils.periodic_structure_analysis import (
    get_template_topology_info,
    calculate_X_X_bonds_and_update_charges,
    adjust_placed_edges_to_fit_periodicity,
    calculate_periodic_distance_and_symmetry,
)
from external_lib.node import Node
from external_lib.edge import Edge
from external_lib.utils.parameter import ShapeMismatchError
from structure_validator import StructureValidator, _get_vdw_radius

class StructureAssembler:
    def __init__(self, output_dir="results/cifs", validator=None, run_symmetry=True,
                 x_site_fusion_tolerance=0.8):
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        self.validator = validator or StructureValidator(strict_mode=False)
        self.run_symmetry = run_symmetry
        self.x_site_fusion_tolerance = float(x_site_fusion_tolerance)
        if self.x_site_fusion_tolerance <= 0:
            raise ValueError("x_site_fusion_tolerance must be positive")

        self.edge_lib_path = os.path.abspath("edges")
        self.edge_lengths = {}
        json_path = os.path.join(self.edge_lib_path, "edge_lengths.json")
        if os.path.exists(json_path):
            with open(json_path, 'r') as f:
                self.edge_lengths = json.load(f)

    def _compute_connection_closure_metrics(self, placed_all, fixed_bonds,
                                            target_bridge_length=1.485,
                                            threshold=0.8):
        """Measure residual node-edge X-site closure errors after PBC-aware bonding."""
        site_by_id = {str(site[0]): site for site in placed_all}
        distances = []

        for bond in fixed_bonds:
            a1, a2 = str(bond[0]), str(bond[1])
            s1 = site_by_id.get(a1)
            s2 = site_by_id.get(a2)
            if s1 is None or s2 is None:
                continue
            if s1[-1] == s2[-1]:
                continue
            if not (str(s1[5]).startswith('X') and str(s2[5]).startswith('X')):
                continue
            distances.append(abs(float(bond[2]) - target_bridge_length))

        if not distances:
            return {
                "n_closure_pairs": 0,
                "closure_rmsd_A": None,
                "closure_median_A": None,
                "closure_p95_A": None,
                "closure_max_A": None,
                "within_0.8A_fraction": None,
                "closure_errors_A": [],
            }

        arr = np.array(distances, dtype=float)
        return {
            "n_closure_pairs": int(len(arr)),
            "closure_rmsd_A": float(np.sqrt(np.mean(arr ** 2))),
            "closure_median_A": float(np.median(arr)),
            "closure_p95_A": float(np.percentile(arr, 95)),
            "closure_max_A": float(np.max(arr)),
            "within_0.8A_fraction": float(np.mean(arr <= threshold)),
            "closure_errors_A": [float(x) for x in arr],
        }

    def _repair_missing_bridge_bonds(self, placed_all, fixed_bonds, sc_unit_cell,
                                     scale=1.0, target_bridge_length=1.485,
                                     closure_threshold=0.8):
        """Conservatively add missing node-edge X bridge bonds after fallback bonding.

        The ToBaCCo-style fallback can under-connect edge X sites when the initial
        graph validation fails. This repair keeps existing bonds and only adds
        one-to-one node↔edge X bonds within the same closure window used by the
        validation trace.
        """
        site_by_id = {str(site[0]): site for site in placed_all}
        node_x = []
        edge_x = []
        for site in placed_all:
            if not str(site[5]).startswith('X'):
                continue
            if site[-1] == 'node':
                node_x.append(site)
            elif site[-1] == 'edge':
                edge_x.append(site)

        used_node_x = set()
        used_edge_x = set()
        existing_bridge_keys = set()
        for bond in fixed_bonds:
            a1, a2 = str(bond[0]), str(bond[1])
            s1 = site_by_id.get(a1)
            s2 = site_by_id.get(a2)
            if s1 is None or s2 is None:
                continue
            if s1[-1] == s2[-1]:
                continue
            if not (str(s1[5]).startswith('X') and str(s2[5]).startswith('X')):
                continue
            existing_bridge_keys.add(tuple(sorted([a1, a2])))
            if s1[-1] == 'node':
                used_node_x.add(a1)
                used_edge_x.add(a2)
            else:
                used_node_x.add(a2)
                used_edge_x.add(a1)

        repair_tol = max(
            (target_bridge_length + closure_threshold) * max(float(scale), 1.0),
            2.75,
        )
        inv_cell = np.linalg.inv(sc_unit_cell)
        candidates = []
        nearest_unmatched = []
        for e_site in edge_x:
            e_id = str(e_site[0])
            if e_id in used_edge_x:
                continue
            e_frac = np.dot(inv_cell, np.array(e_site[1:4], dtype=float))
            nearest_for_edge = None
            for n_site in node_x:
                n_id = str(n_site[0])
                if n_id in used_node_x:
                    continue
                if str(e_site[6]) == str(n_site[6]):
                    continue
                n_frac = np.dot(inv_cell, np.array(n_site[1:4], dtype=float))
                dvec, sym = calculate_periodic_distance_and_symmetry(e_frac, n_frac)
                dist = float(np.linalg.norm(np.dot(sc_unit_cell, dvec)))
                if nearest_for_edge is None or dist < nearest_for_edge["distance"]:
                    nearest_for_edge = {
                        "edge_x": e_id,
                        "node_x": n_id,
                        "distance": dist,
                        "symmetry": sym,
                    }
                if dist <= repair_tol:
                    candidates.append((dist, e_id, n_id, sym))
            if nearest_for_edge is not None:
                nearest_unmatched.append(nearest_for_edge)

        added = []
        candidate_count = len(candidates)
        for dist, e_id, n_id, sym in sorted(candidates, key=lambda item: item[0]):
            key = tuple(sorted([e_id, n_id]))
            if key in existing_bridge_keys:
                continue
            if e_id in used_edge_x or n_id in used_node_x:
                continue
            fixed_bonds.append([e_id, n_id, dist, sym, 'S'])
            used_edge_x.add(e_id)
            used_node_x.add(n_id)
            existing_bridge_keys.add(key)
            added.append({
                "edge_x": e_id,
                "node_x": n_id,
                "distance": dist,
                "symmetry": sym,
            })

        stats = {
            "node_x_count": len(node_x),
            "edge_x_count": len(edge_x),
            "used_node_x_count": len(used_node_x),
            "used_edge_x_count": len(used_edge_x),
            "candidate_pair_count": candidate_count,
            "nearest_unmatched": sorted(nearest_unmatched, key=lambda item: item["distance"])[:12],
            "repair_tolerance_A": repair_tol,
            "added_count": len(added),
        }
        return added, stats

    def _map_final_sites_to_merged_entries(self, structure, merged_placed, tol=0.08):
        """Map final sorted Structure sites back to post-merge assembly entries."""
        if len(structure) == 0:
            return {}

        lattice = structure.lattice.matrix
        inv_lattice = np.linalg.inv(lattice)
        unused = set(range(len(merged_placed)))
        mapped = {}

        for struct_i, site in enumerate(structure):
            specie = str(site.specie)
            site_frac = np.array(site.frac_coords, dtype=float)
            best = None
            for entry_i in list(unused):
                entry = merged_placed[entry_i]
                entry_specie = re.sub('[0-9]', '', str(entry[0]))
                if entry_specie == "X" or entry_specie != specie:
                    continue
                entry_cart = np.array(entry[1:4], dtype=float)
                entry_frac = np.dot(entry_cart, inv_lattice)
                dfrac = site_frac - entry_frac
                dfrac -= np.round(dfrac)
                pbc_dist = float(np.linalg.norm(np.dot(dfrac, lattice)))
                direct_dist = float(np.linalg.norm(site.coords - entry_cart))
                dist = min(pbc_dist, direct_dist)
                if best is None or dist < best[0]:
                    best = (dist, entry_i)
            if best is not None and best[0] <= tol:
                mapped[struct_i] = merged_placed[best[1]]
                unused.remove(best[1])

        return mapped

    def _classify_clash_sources(self, structure, validation_results, merged_placed,
                                merged_bonds, bridge_bonds):
        """Annotate reported clashes with assembly provenance for regression traces."""
        clash = validation_results.get("clash", {}) if isinstance(validation_results, dict) else {}
        reported = clash.get("clashes", []) if isinstance(clash, dict) else []
        if not reported:
            return {"summary": {}, "clashes": []}

        final_to_entry = self._map_final_sites_to_merged_entries(structure, merged_placed)
        bond_keys = {
            tuple(sorted([str(bond[0]), str(bond[1])]))
            for bond in merged_bonds
        }
        bond_graph = nx.Graph()
        for entry in merged_placed:
            bond_graph.add_node(str(entry[0]))
        for bond in merged_bonds:
            a1, a2 = str(bond[0]), str(bond[1])
            if a1 != a2:
                bond_graph.add_edge(a1, a2)

        categories = {}
        label_pairs = {}
        rows = []
        for item in reported:
            i = int(item["i"])
            j = int(item["j"])
            entry_i = final_to_entry.get(i)
            entry_j = final_to_entry.get(j)
            pair_key = None
            path_length = None
            if entry_i is not None and entry_j is not None:
                id_i = str(entry_i[0])
                id_j = str(entry_j[0])
                pair_key = tuple(sorted([id_i, id_j]))
                try:
                    path_length = nx.shortest_path_length(bond_graph, id_i, id_j)
                except (nx.NetworkXNoPath, nx.NodeNotFound):
                    path_length = None

            if entry_i is None or entry_j is None:
                category = "unmapped_after_merge_sites"
            elif pair_key in bond_keys:
                category = "bonded_pair_exclusion_miss"
            elif pair_key in bridge_bonds:
                category = "bridge_bond_exclusion_miss"
            elif str(entry_i[6]) == str(entry_j[6]):
                category = "same_block_near_neighbor" if path_length in (2, 3) else "same_block_nonbonded"
            elif str(entry_i[-1]) == str(entry_j[-1]):
                category = "same_type_interblock"
            else:
                category = "inter_block_nonbonded"

            label_i = str(entry_i[5]) if entry_i is not None else None
            label_j = str(entry_j[5]) if entry_j is not None else None
            label_pair = "-".join(sorted([label_i or "?", label_j or "?"]))
            categories[category] = categories.get(category, 0) + 1
            label_pairs[label_pair] = label_pairs.get(label_pair, 0) + 1

            rows.append({
                "i": i,
                "j": j,
                "elem_i": item.get("elem_i"),
                "elem_j": item.get("elem_j"),
                "distance": item.get("distance"),
                "ratio": item.get("ratio"),
                "category": category,
                "path_length": path_length,
                "same_bbcode": bool(entry_i is not None and entry_j is not None
                                    and str(entry_i[6]) == str(entry_j[6])),
                "same_bbtype": bool(entry_i is not None and entry_j is not None
                                    and str(entry_i[-1]) == str(entry_j[-1])),
                "in_merged_bonds": bool(pair_key in bond_keys) if pair_key is not None else False,
                "in_bridge_bonds": bool(pair_key in bridge_bonds) if pair_key is not None else False,
                "within_merge_sites_tol": bool(item.get("distance", 999.0) <= 0.2),
                "near_merge_sites_tol": bool(item.get("distance", 999.0) <= 0.25),
                "atom_i": None if entry_i is None else {
                    "id": str(entry_i[0]),
                    "label": str(entry_i[5]),
                    "bbcode": str(entry_i[6]),
                    "bbtype": str(entry_i[-1]),
                },
                "atom_j": None if entry_j is None else {
                    "id": str(entry_j[0]),
                    "label": str(entry_j[5]),
                    "bbcode": str(entry_j[6]),
                    "bbtype": str(entry_j[-1]),
                },
            })

        return {
            "summary": {
                "category_counts": dict(sorted(categories.items())),
                "top_label_pairs": sorted(
                    label_pairs.items(), key=lambda kv: (-kv[1], kv[0])
                )[:12],
                "mapped_clash_atoms": sum(
                    1 for row in rows if row["atom_i"] is not None and row["atom_j"] is not None
                ),
                "reported_clashes": len(rows),
                "total_clashes": clash.get("n_clashes"),
            },
            "clashes": rows,
        }

    def _axis_rotation_matrix(self, axis, angle_rad):
        axis = np.array(axis, dtype=float)
        norm = np.linalg.norm(axis)
        if norm < 1e-12:
            return np.eye(3)
        axis = axis / norm
        x, y, z = axis
        c = np.cos(angle_rad)
        s = np.sin(angle_rad)
        C = 1.0 - c
        return np.array([
            [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
        ])

    def _relieve_edge_axial_clashes(self, placed_all, fixed_bonds, sc_unit_cell,
                                    clash_tolerance=0.6):
        """Rotate 2-connected edge blocks around their X-X axis to reduce contacts."""
        edge_blocks = {}
        for idx, site in enumerate(placed_all):
            if site[-1] == 'edge':
                edge_blocks.setdefault(str(site[6]), []).append(idx)
        if not edge_blocks:
            return []

        atom_to_idx = {str(site[0]): idx for idx, site in enumerate(placed_all)}
        bond_graph = nx.Graph()
        for site in placed_all:
            bond_graph.add_node(str(site[0]))
        for bond in fixed_bonds:
            a1, a2 = str(bond[0]), str(bond[1])
            if a1 in atom_to_idx and a2 in atom_to_idx and a1 != a2:
                bond_graph.add_edge(a1, a2)

        exclusion_keys = set()
        for atom_id in bond_graph:
            lengths = nx.single_source_shortest_path_length(bond_graph, atom_id, cutoff=2)
            for other_id, depth in lengths.items():
                if depth > 0:
                    i = atom_to_idx.get(atom_id)
                    j = atom_to_idx.get(other_id)
                    if i is not None and j is not None and i != j:
                        exclusion_keys.add(tuple(sorted((i, j))))

        inv_cell = np.linalg.inv(sc_unit_cell)
        species = [re.sub('[0-9]', '', str(site[0])) for site in placed_all]
        radii = [_get_vdw_radius(sp) for sp in species]

        def _coords():
            return np.array([list(map(float, site[1:4])) for site in placed_all], dtype=float)

        def _pbc_distance(coords, i, j):
            fi = np.dot(inv_cell, coords[i])
            fj = np.dot(inv_cell, coords[j])
            dfrac = fi - fj
            dfrac -= np.round(dfrac)
            return float(np.linalg.norm(np.dot(sc_unit_cell, dfrac)))

        def _score_block(coords, block_indices):
            score = 0.0
            count = 0
            block_set = set(block_indices)
            for i in block_indices:
                if species[i] == "X":
                    continue
                for j in range(len(coords)):
                    if j == i or j in block_set or species[j] == "X":
                        continue
                    key = tuple(sorted((i, j)))
                    if key in exclusion_keys:
                        continue
                    vdw_sum = radii[i] + radii[j]
                    if vdw_sum <= 0:
                        continue
                    ratio = _pbc_distance(coords, i, j) / vdw_sum
                    if ratio < clash_tolerance:
                        count += 1
                        score += 1.0 + (clash_tolerance - ratio) ** 2 * 100.0
                    elif ratio < 0.72:
                        score += (0.72 - ratio) ** 2
            return score, count

        coords = _coords()
        trace = []
        angles = [np.deg2rad(deg) for deg in range(0, 360, 30)]
        for bbcode, block_indices in sorted(edge_blocks.items(), key=lambda item: int(item[0])):
            anchors = [
                idx for idx in block_indices
                if str(placed_all[idx][5]).startswith('X')
            ]
            if len(anchors) != 2:
                continue
            p0 = coords[anchors[0]]
            p1 = coords[anchors[1]]
            axis = p1 - p0
            if np.linalg.norm(axis) < 1e-8:
                continue
            base_score, base_count = _score_block(coords, block_indices)
            best = (base_score, base_count, 0.0, coords[block_indices].copy())
            for angle in angles[1:]:
                rot = self._axis_rotation_matrix(axis, angle)
                trial_coords = coords.copy()
                shifted = coords[block_indices] - p0
                trial_coords[block_indices] = np.dot(shifted, rot.T) + p0
                score, count = _score_block(trial_coords, block_indices)
                if (count, score) < (best[1], best[0]):
                    best = (score, count, angle, trial_coords[block_indices].copy())
            if best[2] != 0.0 and (best[1], best[0]) < (base_count, base_score):
                coords[block_indices] = best[3]
                angle_deg = float(np.rad2deg(best[2]))
                trace.append({
                    "edge_bbcode": bbcode,
                    "angle_deg": angle_deg,
                    "score_before": float(base_score),
                    "score_after": float(best[0]),
                    "clashes_before": int(base_count),
                    "clashes_after": int(best[1]),
                })

        if trace:
            for i, site in enumerate(placed_all):
                site[1:4] = [float(x) for x in coords[i]]
            print(f"  [EdgeRelax] Rotated {len(trace)} edge blocks to reduce axial clashes")
        return trace

    def find_best_edge(self, target_length):
        if not self.edge_lengths: return None, None
        best_name = None
        min_diff = float('inf')
        for name, length in self.edge_lengths.items():
            diff = abs(length - target_length)
            if diff < min_diff:
                min_diff = diff
                best_name = name
        return best_name, min_diff

    def get_node_radius(self, cif_name):
        if not hasattr(self, '_radius_cache'):
            self._radius_cache = {}
        if cif_name not in self._radius_cache:
            try:
                temp_node = Node(os.path.join('nodes', f"{cif_name}.cif"))
                x_coords = temp_node.calculate_X_cart_vectors(False)
                if len(x_coords) > 0:
                    self._radius_cache[cif_name] = np.mean([np.linalg.norm(v) for v in x_coords])
                else:
                    self._radius_cache[cif_name] = 3.5
            except Exception:
                self._radius_cache[cif_name] = 3.5
        return self._radius_cache[cif_name]

    def _get_edge_span(self, cif_name):
        if not hasattr(self, '_edge_span_cache'):
            self._edge_span_cache = {}
        if cif_name in self._edge_span_cache:
            return self._edge_span_cache[cif_name]
        try:
            temp_edge = Edge(os.path.join('edges', f"{cif_name}.cif"))
            unit_cell, fcoords = temp_edge.calculate_unitcell_fcoords()
            if len(fcoords) < 2:
                self._edge_span_cache[cif_name] = 0.0
                return 0.0
            cart = [np.dot(unit_cell, f[1]) for f in fcoords]
            max_d = max(np.linalg.norm(cart[i] - cart[j])
                        for i in range(len(cart)) for j in range(i + 1, len(cart)))
            self._edge_span_cache[cif_name] = max_d
            return max_d
        except Exception:
            self._edge_span_cache[cif_name] = 0.0
            return 0.0

    def _compute_lattice_scale(self, graph_data, structure_map, forced_edge=None):
        lattice_mat = np.array(graph_data['lattice'])
        nodes_info = graph_data['nodes']
        edges_data = graph_data['edges']
        default_Le = np.mean(list(self.edge_lengths.values())) if self.edge_lengths else 5.0
        max_scale = 1.0
        for edge in edges_data:
            u, v = str(edge['source']), str(edge['target'])
            R_u = self.get_node_radius(structure_map[u])
            R_v = self.get_node_radius(structure_map[v])

            d_frac = np.array(nodes_info[edge['target']]['frac']) - np.array(nodes_info[edge['source']]['frac'])
            d_frac = d_frac - np.round(d_frac)
            actual_dist = np.linalg.norm(np.dot(lattice_mat, d_frac))
            actual_Le = actual_dist - R_u - R_v
            if actual_Le < 0.1:
                actual_Le = default_Le

            if forced_edge:
                edge_cif = forced_edge
            else:
                edge_cif, _ = self.find_best_edge(actual_Le)
                if not edge_cif:
                    continue

            L_edge = self._get_edge_span(edge_cif)
            if L_edge <= 0:
                continue

            required = R_u + L_edge + R_v
            if actual_dist > 0 and actual_dist < required:
                max_scale = max(max_scale, required / actual_dist)

        if max_scale > 1.0:
            max_scale *= 1.05
        return max_scale

    def assemble_and_save(self, filename, structure_map, graph_data, forced_edge=None,
                          return_metadata=False):
        lattice_mat = np.array(graph_data['lattice'])
        nodes_info = graph_data['nodes']
        edges_data = graph_data['edges']

        # Expand the unit cell when needed so the chosen linker can fit geometrically.
        scale = self._compute_lattice_scale(graph_data, structure_map, forced_edge)
        if abs(scale - 1.0) > 1e-6:
            print(f"  [Lattice] Scaling unit cell by {scale:.4f}x to accommodate ligand sizes")
            lattice_mat = lattice_mat * scale
        else:
            print(f"  [Lattice] No scaling needed (scale={scale:.4f})")

        num_vertices = len(nodes_info)
        num_edges = len(edges_data)

        print("  [Topology] Building NetworkX graph and geometry definitions...")

        TG = nx.MultiGraph()
        for i, node in enumerate(nodes_info):
            node_id = str(node.get('id', i))
            cifname = f"{structure_map[node_id]}.cif"
            TG.add_node(node_id, cifname=cifname)

        ea_dict = {str(node.get('id', i)): {} for i, node in enumerate(nodes_info)}
        default_Le = np.mean(list(self.edge_lengths.values())) if self.edge_lengths else 5.0

        for i, edge in enumerate(edges_data):
            u, v = str(edge['source']), str(edge['target'])
            R_u = self.get_node_radius(structure_map[u])
            R_v = self.get_node_radius(structure_map[v])

            d_frac = np.array(nodes_info[edge['target']]['frac']) - np.array(nodes_info[edge['source']]['frac'])
            actual_distance = np.linalg.norm(np.dot(lattice_mat, d_frac))
            actual_Le = actual_distance - R_u - R_v
            if actual_Le < 0.1: actual_Le = default_Le

            if forced_edge:
                edge_cif = forced_edge
                print(f"  [Assembler] Using forced edge: {edge_cif}")
            else:
                edge_cif, _ = self.find_best_edge(actual_Le)
                if not edge_cif: edge_cif = list(self.edge_lengths.keys())[0]
            edge_index = i + 1
            TG.add_edge(u, v, index=edge_index, pd=(u, v), cifname=f"{edge_cif}.cif")

            ea_dict[u][edge_index] = [None, R_u]
            ea_dict[v][edge_index] = [None, R_v]

        prealign_trace = []
        for i, node in enumerate(nodes_info):
            node_id = str(node.get('id', i))
            node_cif = structure_map[node_id]
            try:
                node_x_count = len(Node(os.path.join('nodes', f"{node_cif}.cif")).calculate_X_cart_vectors(False))
            except Exception:
                node_x_count = None
            anchor_count = len(ea_dict.get(node_id, {}))
            prealign_trace.append({
                "node": node_id,
                "cif": node_cif,
                "node_x_count": node_x_count,
                "topology_anchor_count": anchor_count,
            })
            print(
                "  [PreAlign] "
                f"node={node_id} cif={node_cif} "
                f"node_x_count={node_x_count} topology_anchor_count={anchor_count}"
            )

        l = Lattice(lattice_mat.T)
        sc_abc = (l.a, l.b, l.c)
        sc_Roman = (l.alpha, l.beta, l.gamma)
        start = nodes_info[0]['frac']

        alpha_matrix = np.array(graph_data.get('alpha_matrix', np.zeros((num_edges - num_vertices + 1, 3))))
        Bstar = np.array(graph_data.get('Bstar', np.eye(3)))
        Bstar_inv = np.linalg.pinv(Bstar)
        sc_covar = np.zeros((num_vertices - 1, 3))

        print("  [Alignment] Unrolling periodic topology using exact mathematical scaling...")
        try:
            (sc_unit_cell, placed_nodes, node_bonds,
             placed_edges, edge_bonds, modified, removed) = get_template_topology_info(
                alpha_matrix, num_edges, num_vertices, sc_covar, Bstar_inv,
                sc_abc, sc_Roman, start, TG, "template", None, ea_dict
            )
        except ShapeMismatchError as e:
            print(f"  [ShapeMismatch] {e}")
            print(f"  [Assembly] Skipping this topology-ligand combination due to shape mismatch.")
            raise
        except Exception as e:
            print(f"  [Error] Original topology solver failed: {e}")
            raise e

        placed_edges = adjust_placed_edges_to_fit_periodicity(placed_edges, placed_nodes, sc_unit_cell)

        placed_all = [list(n) + ['node'] for n in placed_nodes] + [list(e) + ['edge'] for e in placed_edges]
        bonds_all = node_bonds + edge_bonds

        blocks = {}
        for n in placed_all:
            if n[6] not in blocks: blocks[n[6]] = []
            blocks[n[6]].append(n[0])
        existing_bonds = set((b[0], b[1]) for b in bonds_all) | set((b[1], b[0]) for b in bonds_all)
        for bbind, atoms in blocks.items():
            if len(atoms) > 1:
                for leaf in atoms[1:]:
                    if (atoms[0], leaf) not in existing_bonds:
                        bonds_all.append([atoms[0], leaf, 1.0, '.', 'S'])

        scaled_max_length = 5.0 * scale
        scaled_bond_tol = 2.5 * scale
        fixed_bonds, fc_placed_all, _ = calculate_X_X_bonds_and_update_charges(
            placed_all, bonds_all, sc_unit_cell, max_length=scaled_max_length,
            nconnections=[[n[0], 1] for n in placed_all if "X" in str(n[5])],
            num_possible_XX_bonds=num_edges * 2,
            bond_tol=scaled_bond_tol
        )

        bridge_repair_trace, bridge_repair_stats = self._repair_missing_bridge_bonds(
            fc_placed_all, fixed_bonds, sc_unit_cell, scale=scale
        )
        if bridge_repair_trace:
            print(f"  [BridgeRepair] Added {len(bridge_repair_trace)} missing node↔edge X bonds")

        edge_rotation_trace = self._relieve_edge_axial_clashes(
            fc_placed_all, fixed_bonds, sc_unit_cell,
            clash_tolerance=self.validator.clash_tolerance,
        )

        # Merge overlapping node-side and edge-side bridge atoms into a single site.
        fragment_of = {str(e[0]): e[-1] for e in fc_placed_all}
        closure_metrics = self._compute_connection_closure_metrics(fc_placed_all, fixed_bonds)

        merge_map = {}
        matched_edge_ids = set()
        for ni in fc_placed_all:
            if ni[-1] != 'node' or not str(ni[5]).startswith('X'):
                continue
            n_pos = np.array(ni[1:4], dtype=float)
            for ej in fc_placed_all:
                if ej[-1] != 'edge' or not str(ej[5]).startswith('X'):
                    continue
                if ej[0] in matched_edge_ids:
                    continue
                e_pos = np.array(ej[1:4], dtype=float)
                if np.linalg.norm(n_pos - e_pos) < self.x_site_fusion_tolerance:
                    merge_map[ni[0]] = ej[0]
                    matched_edge_ids.add(ej[0])
                    break

        topology_linked_x_pairs = set()
        placed_by_id = {str(site[0]): site for site in fc_placed_all}
        for bond in fixed_bonds:
            a1, a2 = str(bond[0]), str(bond[1])
            s1, s2 = placed_by_id.get(a1), placed_by_id.get(a2)
            if (s1 is not None and s2 is not None and s1[-1] != s2[-1]
                    and str(s1[5]).startswith('X') and str(s2[5]).startswith('X')):
                node_id, edge_id = (a1, a2) if s1[-1] == 'node' else (a2, a1)
                topology_linked_x_pairs.add((node_id, edge_id))
        selected_x_pairs = set(merge_map.items())
        fusion_audit = {
            "tolerance_A": self.x_site_fusion_tolerance,
            "topology_linked_pairs": len(topology_linked_x_pairs),
            "merged_pairs": len(selected_x_pairs),
            "topology_linked_pairs_merged": len(topology_linked_x_pairs & selected_x_pairs),
            "unlinked_pairs_merged": len(selected_x_pairs - topology_linked_x_pairs),
            "topology_linked_pair_ids": sorted(topology_linked_x_pairs),
            "merged_pair_ids": sorted(selected_x_pairs),
        }

        removed_node_ids = set(merge_map.keys())
        merged_placed = [site for site in fc_placed_all if site[0] not in removed_node_ids]

        bridge_bonds = set()
        merged_bonds_raw = []
        for bond in fixed_bonds:
            a1, a2 = str(bond[0]), str(bond[1])
            if a1 in merge_map:
                a1 = merge_map[a1]
            if a2 in merge_map:
                a2 = merge_map[a2]
            if a1 == a2:
                continue
            dist_f = float(bond[2])
            if dist_f < 0.8:
                continue
            f1 = fragment_of.get(str(bond[0]))
            f2 = fragment_of.get(str(bond[1]))
            if f1 is not None and f2 is not None and f1 != f2:
                bridge_bonds.add(tuple(sorted([str(a1), str(a2)])))
            merged_bonds_raw.append([a1, a2, bond[2], bond[3], bond[4]])

        seen = set()
        merged_bonds = []
        for b in merged_bonds_raw:
            key = tuple(sorted([str(b[0]), str(b[1])]))
            if key not in seen:
                seen.add(key)
                merged_bonds.append(b)

        if merge_map:
            print(f"  [Merge] Merged {len(merge_map)} node-side X atoms into edge-side bridge C atoms")
        print(f"  [Bridge] Identified {len(bridge_bonds)} Node↔Edge bridge bonds for 1.485 Å overwrite")

        species = [re.sub('[0-9]', '', str(site[0])) for site in merged_placed]
        coords = [list(map(float, site[1:4])) for site in merged_placed]

        final_struct = Structure(Lattice(sc_unit_cell.T), species, coords, coords_are_cartesian=True)
        atom_trace = {
            "merged_placed_sites": len(merged_placed),
            "x_sites_before_remove": sum(1 for sp in species if sp == "X"),
            "sites_before_remove_x": len(final_struct),
        }
        final_struct.remove_species(["X"])
        atom_trace["sites_after_remove_x"] = len(final_struct)

        atom_trace["sites_before_merge_sites"] = len(final_struct)
        final_struct.merge_sites(tol=0.2, mode="delete")
        atom_trace["sites_after_merge_sites"] = len(final_struct)

        detected_sg = "P1"
        atom_trace["sites_before_symmetry_refine"] = len(final_struct)
        if self.run_symmetry:
            try:
                sga = SpacegroupAnalyzer(final_struct, symprec=0.3)
                sg_symbol = sga.get_space_group_symbol()

                if sg_symbol:
                    detected_sg = sg_symbol
                print(f"  [Symmetry] Symmetrized structure Spacegroup: {detected_sg}")

                if detected_sg != "P1":
                    refined = sga.get_refined_structure()
                    if refined:
                        final_struct = refined
            except TypeError:
                print(
                    f"  [Symmetry] spglib symmetry detection aborted (likely due to redundant edge artifacts). Fallback to P1.")
            except Exception as e:
                print(f"  [Symmetry] Spacegroup analysis bypassed: {e}. Fallback to P1.")
        else:
            print("  [Symmetry] Spacegroup analysis skipped for reconstruction audit.")
        atom_trace["sites_after_symmetry_refine"] = len(final_struct)

        final_struct = final_struct.get_sorted_structure()
        atom_trace["sites_after_sort"] = len(final_struct)
        print(
            "  [AtomTrace] "
            + " ".join(f"{key}={value}" for key, value in atom_trace.items())
        )

        validation_results = self.validator.validate(
            structure=final_struct,
            topology_graph=TG,
            placed_all=fc_placed_all,
            bonds_all=fixed_bonds,
            lattice_matrix=sc_unit_cell,
            clash_placed_all=merged_placed,
            clash_bonds_all=merged_bonds,
        )
        clash_source_trace = self._classify_clash_sources(
            final_struct, validation_results, merged_placed, merged_bonds, bridge_bonds
        )

        if not validation_results['passed']:
            print(f"  [Validation] WARNING: Structure validation failed!")
            for check_name in ('clash', 'topology', 'symmetry'):
                cr = validation_results.get(check_name, {})
                if not cr.get('skipped'):
                    status = "PASS" if cr.get('passed', False) else "FAIL"
                    print(f"    [{status}] {check_name}: {cr.get('message', cr.get('detected_sg', '?'))}")
        else:
            print(f"  [Validation] All checks passed.")
            if validation_results.get('symmetry', {}).get('detected_sg'):
                print(f"    [Symmetry] {validation_results['symmetry']['message']}")

        final_struct.to(filename=filename)

        # Rewrite atom labels to a 1-indexed CIF-friendly convention.
        sorted_sites = final_struct.sites
        atom_id_to_label = {}
        for entry in merged_placed:
            atom_id = str(entry[0])
            ec = np.array(entry[1:4], dtype=float)
            for j, site in enumerate(sorted_sites):
                if np.linalg.norm(ec - site.coords) < 0.01:
                    atom_id_to_label[atom_id] = f"{site.specie.symbol}{j}"
                    break

        def _label_to_1indexed(m):
            return f"{m.group(1)}{int(m.group(2)) + 1}"

        with open(filename, 'r') as f:
            content = f.read()
        lines = content.split('\n')
        new_lines = []
        for line in lines:
            if line and line[0] in ('_', 'l', 'd', '#'):
                new_lines.append(line)
            else:
                new_lines.append(re.sub(r'\b([A-Z][a-z]?)(\d+)\b', _label_to_1indexed, line))
        with open(filename, 'w') as f:
            f.write('\n'.join(new_lines))

        with open(filename, 'a') as f:
            f.write("\nloop_\n _geom_bond_atom_site_label_1\n")
            f.write(" _geom_bond_atom_site_label_2\n")
            f.write(" _geom_bond_distance\n")
            f.write(" _geom_bond_site_symmetry_2\n")
            f.write(" _ccdc_geom_bond_type\n")
            for bond in merged_bonds:
                a1, a2, dist, symmetry, bond_type = bond
                dist_f = float(dist)
                if abs(dist_f - 1.0) < 0.001:
                    continue
                label1 = atom_id_to_label.get(str(a1))
                label2 = atom_id_to_label.get(str(a2))
                if label1 is None or label2 is None:
                    continue
                label1 = re.sub(r'\b([A-Z][a-z]?)(\d+)\b', _label_to_1indexed, label1)
                label2 = re.sub(r'\b([A-Z][a-z]?)(\d+)\b', _label_to_1indexed, label2)
                bond_key = tuple(sorted([str(a1), str(a2)]))
                if bond_key in bridge_bonds:
                    dist = 1.485
                    bond_type = 'S'
                f.write(f"{label1:>8s} {label2:>8s} {float(dist):>6.3f}   {symmetry:6s} {bond_type}\n")

        print(f"  [Assembly] Assembly Completed. Saved to: {filename}")

        if return_metadata:
            return {
                "cif_path": filename,
                "scale_factor": float(scale),
                "lattice_scaled": abs(scale - 1.0) > 1e-6,
                "merge_count": len(merge_map),
                "bridge_bond_count": len(bridge_bonds),
                "bridge_repair_count": len(bridge_repair_trace),
                "bridge_repair_trace": bridge_repair_trace,
                "bridge_repair_stats": bridge_repair_stats,
                "edge_rotation_trace": edge_rotation_trace,
                "connection_closure": closure_metrics,
                "x_site_fusion": fusion_audit,
                "detected_sg": detected_sg,
                "n_atoms": len(final_struct),
                "prealign_trace": prealign_trace,
                "atom_trace": atom_trace,
                "clash_source_trace": clash_source_trace,
                "validation": validation_results,
                "validation_passed": bool(validation_results.get('passed', False)),
            }
