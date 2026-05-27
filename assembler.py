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
    adjust_placed_edges_to_fit_periodicity
)
from external_lib.node import Node
from external_lib.edge import Edge
from external_lib.utils.parameter import ShapeMismatchError
from structure_validator import StructureValidator


class StructureAssembler:
    def __init__(self, output_dir="results/cifs", validator=None):
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        self.validator = validator or StructureValidator(strict_mode=False)

        self.edge_lib_path = os.path.abspath("edges")
        self.edge_lengths = {}
        json_path = os.path.join(self.edge_lib_path, "edge_lengths.json")
        if os.path.exists(json_path):
            with open(json_path, 'r') as f:
                self.edge_lengths = json.load(f)

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

    def assemble_and_save(self, filename, structure_map, graph_data, forced_edge=None):
        lattice_mat = np.array(graph_data['lattice'])
        nodes_info = graph_data['nodes']
        edges_data = graph_data['edges']

        # ── 晶胞自适应缩放：防止配体间位阻冲突 ──
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
            # ----------------------
            edge_index = i + 1
            TG.add_edge(u, v, index=edge_index, pd=(u, v), cifname=f"{edge_cif}.cif")

            ea_dict[u][edge_index] = [None, R_u]
            ea_dict[v][edge_index] = [None, R_v]

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

        # 6. 生成结构与空间群恢复
        # ── 连接点原子融合（Merge）：代替 Phase 2A 一刀切删除 ──
        # 配体端 X2/X5（type C）和节点端 X11-X14（type C）在连接点处重叠（~0.6 Å），
        # 此处将配对融合为单个真实 C 原子，保留配体侧坐标并移除节点侧副本，
        # 节点侧键通过 merge_map 重定向到保留的原子。

        # ── 片段归属表：用于识别 Node↔Edge 桥接键 ──
        fragment_of = {str(e[0]): e[-1] for e in fc_placed_all}

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
                if np.linalg.norm(n_pos - e_pos) < 0.8:
                    merge_map[ni[0]] = ej[0]
                    matched_edge_ids.add(ej[0])
                    break

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
                continue  # 跳过融合后自环
            dist_f = float(bond[2])
            if dist_f < 0.8:
                continue  # 跳过融合前的 X-X 短键
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
        final_struct.remove_species(["X"])

        # --- 【防弹核心补丁：自愈合清理】 ---
        # 容差必须设置为 0.2 左右，1.0 会导致正常的 C-C 键和 C-H 键原子被误删
        final_struct.merge_sites(tol=0.2, mode="delete")

        detected_sg = "P1"
        try:
            # 保持较大的 symprec 以应对微小形变
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

        final_struct = final_struct.get_sorted_structure()

        # ── 结构自校验 ──
        validation_results = self.validator.validate(
            structure=final_struct,
            topology_graph=TG,
            placed_all=fc_placed_all,
            bonds_all=fixed_bonds,
            lattice_matrix=sc_unit_cell,
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

        # ── CIF 后处理：1-indexed 标签重命名 + 追加 _geom_bond 表 ──
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
                    continue  # 跳过未更新的占位键（距离 1.0 Å 为 block 连接代码遗留的占位符）
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