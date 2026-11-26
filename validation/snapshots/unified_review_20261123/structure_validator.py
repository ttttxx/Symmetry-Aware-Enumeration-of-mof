"""
StructureValidator — MOF 3D 结构自校验引擎

提供三层硬性校验:
  1. 空间位阻/碰撞检测 (Clash Detection)
  2. 拓扑还原校验 (Topology Isomorphism Check)
  3. 对称性守恒校验 (Symmetry Retention Check)
"""
import logging
from collections import defaultdict

import networkx as nx
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

logger = logging.getLogger(__name__)


# ── 常见元素范德华半径 (Å) ────────────────────────────────────────────
VDW_RADII = {
    'H': 1.20, 'He': 1.40, 'Li': 1.82, 'Be': 1.53, 'B': 1.92,
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'Ne': 1.54,
    'Na': 2.27, 'Mg': 1.73, 'Al': 1.84, 'Si': 2.10, 'P': 1.80,
    'S': 1.80, 'Cl': 1.75, 'Ar': 1.88, 'K': 2.75, 'Ca': 2.31,
    'Sc': 2.11, 'Ti': 2.11, 'V': 2.11, 'Cr': 2.11, 'Mn': 2.11,
    'Fe': 2.11, 'Co': 2.11, 'Ni': 2.11, 'Cu': 2.11, 'Zn': 2.11,
    'Ga': 1.87, 'Ge': 1.96, 'As': 1.85, 'Se': 1.90, 'Br': 1.85,
    'Kr': 2.02, 'Rb': 3.03, 'Sr': 2.49, 'Y': 2.31, 'Zr': 2.31,
    'Nb': 2.31, 'Mo': 2.31, 'Tc': 2.31, 'Ru': 2.31, 'Rh': 2.31,
    'Pd': 2.31, 'Ag': 2.31, 'Cd': 2.31, 'In': 1.93, 'Sn': 2.17,
    'Sb': 2.06, 'Te': 2.06, 'I': 1.98, 'Xe': 2.16, 'Cs': 3.43,
    'Ba': 2.68, 'La': 2.43, 'Ce': 2.43, 'Pr': 2.43, 'Nd': 2.43,
    'Pm': 2.43, 'Sm': 2.43, 'Eu': 2.43, 'Gd': 2.43, 'Tb': 2.43,
    'Dy': 2.43, 'Ho': 2.43, 'Er': 2.43, 'Tm': 2.43, 'Yb': 2.43,
    'Lu': 2.43, 'Hf': 2.31, 'Ta': 2.31, 'W': 2.31, 'Re': 2.31,
    'Os': 2.31, 'Ir': 2.31, 'Pt': 2.31, 'Au': 2.31, 'Hg': 2.31,
    'Tl': 1.96, 'Pb': 2.02, 'Bi': 2.07, 'Po': 1.97, 'At': 2.02,
    'Rn': 2.20, 'Fr': 3.48, 'Ra': 2.83, 'Ac': 2.43, 'Th': 2.43,
    'Pa': 2.43, 'U': 2.43, 'Np': 2.43, 'Pu': 2.43, 'Am': 2.43,
    'Cm': 2.43, 'Bk': 2.43, 'Cf': 2.43, 'Es': 2.43, 'Fm': 2.43,
    'Md': 2.43, 'No': 2.43, 'Lr': 2.43,
}

# ── 常见元素共价半径 (Å)，用于拓扑还原中的成键检测 ──────────────────
COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07,
    'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76,
    'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39,
    'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
    'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20,
    'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75,
    'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42,
    'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
    'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40, 'Cs': 2.44,
    'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01,
    'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94,
    'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
    'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62, 'Re': 1.51,
    'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
    'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 'At': 1.50,
    'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
    'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87, 'Am': 1.80,
    'Cm': 1.69,
}


def _get_vdw_radius(element):
    return VDW_RADII.get(element, 1.80)


def _get_covalent_radius(element):
    return COVALENT_RADII.get(element, 0.76)


class StructureValidator:
    """MOF 结构自校验器。

    Args:
        clash_tolerance: 碰撞阈值系数 (默认 0.6, 即 d < 0.6 * (vdW_i + vdW_j) 视为碰撞)
        bond_tolerance:  成键检测容差 (默认 0.5 Å, 加到共价半径和上)
        symprec:         对称性检测精度 (默认 0.3)
        strict_mode:     严格模式 — 任何校验失败都标记 passed=False
    """

    def __init__(self, clash_tolerance=0.6, bond_tolerance=0.5,
                 symprec=0.3, strict_mode=True, check_symmetry=True):
        self.clash_tolerance = clash_tolerance
        self.bond_tolerance = bond_tolerance
        self.symprec = symprec
        self.strict_mode = strict_mode
        self.check_symmetry = check_symmetry

    # ================================================================
    #  公开接口
    # ================================================================

    def validate(self, structure, topology_graph=None, placed_all=None,
                 bonds_all=None, lattice_matrix=None, clash_placed_all=None,
                 clash_bonds_all=None, clash_bond_graph_depth=2):
        """执行三层校验。

        Args:
            structure:      pymatgen Structure 对象 (应已去除 X 原子)
            topology_graph: NetworkX Graph，原始拓扑图
            placed_all:     放置后的所有原子信息 (含 bbcode/bbtype)，用于拓扑校验
            bonds_all:      所有化学键列表，用于拓扑校验
            lattice_matrix: 晶格矩阵 (3x3)，用于拓扑校验
            clash_placed_all: 可选的碰撞排除原子表；默认使用 placed_all
            clash_bonds_all:  可选的碰撞排除键表；默认使用 bonds_all
            clash_bond_graph_depth: 排除几何 clash 的键图距离；2 表示排除直接键和 1-3 邻接

        Returns:
            dict: 包含 'passed' (bool) 和各检查详细结果
        """
        results = {}

        # ── 碰撞检测的键合排除集 ──
        # 从 bonds_all 和 placed_all 中提取键合原子对的索引，
        # 通过坐标匹配关联到 structure 中的位点。
        clash_placed = clash_placed_all if clash_placed_all is not None else placed_all
        clash_bonds = clash_bonds_all if clash_bonds_all is not None else bonds_all
        bond_exclusions = self._build_clash_exclusions(
            structure, clash_placed, clash_bonds, max_graph_depth=clash_bond_graph_depth
        )

        results['clash'] = self._check_clashes(structure, bond_exclusions)

        # ── 2. 拓扑还原校验 ──
        if topology_graph is not None and placed_all is not None:
            results['topology'] = self._check_topology(
                structure, topology_graph, placed_all, bonds_all, lattice_matrix
            )
        else:
            results['topology'] = {
                'passed': True,
                'skipped': True,
                'message': '缺少拓扑数据，跳过校验'
            }

        # ── 3. 对称性守恒校验 ──
        if self.check_symmetry:
            results['symmetry'] = self._check_symmetry(
                structure, topology_graph
            )
        else:
            results['symmetry'] = {
                'passed': True,
                'skipped': True,
                'detected_sg': 'not evaluated',
                'message': '对称性分析已按本次审计设置跳过',
            }

        # ── 综合判定 ──
        if self.strict_mode:
            results['passed'] = all(
                v.get('passed', False) for v in results.values()
            )
        else:
            # 非严格模式下仅 clash 和 topology 的失败算失败
            fails = []
            for key in ('clash', 'topology'):
                v = results.get(key, {})
                if v.get('passed') is False and not v.get('skipped'):
                    fails.append(key)
            results['passed'] = len(fails) == 0

        return results

    def _build_clash_exclusions(self, structure, placed_all, bonds_all, max_graph_depth=2):
        """Map placed bond graph neighborhoods to final Structure indices."""
        if bonds_all is None or placed_all is None or len(placed_all) == 0:
            return None

        lattice = structure.lattice.matrix
        inv_lattice = np.linalg.inv(lattice)
        placed_to_struct = {}
        atom_to_placed = {}
        for i, entry in enumerate(placed_all):
            atom_to_placed[str(entry[0])] = i
            placed_specie = ''.join(ch for ch in str(entry[0]) if not ch.isdigit())
            if placed_specie == "X":
                continue
            pc = np.array(list(map(float, entry[1:4])), dtype=float)
            pc_frac = np.dot(pc, inv_lattice)
            for j, site in enumerate(structure):
                if str(site.specie) != placed_specie:
                    continue
                dfrac = np.array(site.frac_coords, dtype=float) - pc_frac
                dfrac -= np.round(dfrac)
                pbc_dist = np.linalg.norm(np.dot(dfrac, lattice))
                direct_dist = np.linalg.norm(pc - np.array(site.coords))
                if min(pbc_dist, direct_dist) < 0.08:
                    placed_to_struct[i] = j
                    break

        atom_graph = nx.Graph()
        for atom_id in atom_to_placed:
            atom_graph.add_node(atom_id)
        for bond in bonds_all:
            a1, a2 = str(bond[0]), str(bond[1])
            if a1 in atom_to_placed and a2 in atom_to_placed and a1 != a2:
                atom_graph.add_edge(a1, a2)

        bond_exclusions = set()
        if max_graph_depth is None or max_graph_depth < 1:
            max_graph_depth = 1
        for atom_id in atom_graph:
            pi = atom_to_placed.get(atom_id)
            si = placed_to_struct.get(pi)
            if si is None:
                continue
            lengths = nx.single_source_shortest_path_length(
                atom_graph, atom_id, cutoff=max_graph_depth
            )
            for other_id, depth in lengths.items():
                if depth == 0:
                    continue
                pj = atom_to_placed.get(other_id)
                sj = placed_to_struct.get(pj)
                if sj is None or si == sj:
                    continue
                bond_exclusions.add(tuple(sorted((si, sj))))

        return bond_exclusions

    # ================================================================
    #  碰撞检测
    # ================================================================

    def _check_clashes(self, structure, exclude_bonds=None):
        """遍历所有非键合原子对，检测空间位阻。

        判断标准: d < clash_tolerance × (vdW_i + vdW_j)

        Args:
            structure:     pymatgen Structure 对象
            exclude_bonds: 可选的 set，包含应排除的 (i, j) 键合原子对索引
        """
        clashes = []
        min_ratio = float('inf')
        n_sites = len(structure)

        if n_sites < 2:
            return {
                'passed': True,
                'clashes': [],
                'min_distance': None,
                'message': '原子数不足 2，跳过碰撞检测'
            }

        # 预计算每个位点的 vdW 半径
        radii = [_get_vdw_radius(str(site.specie)) for site in structure]

        for i in range(n_sites):
            for j in range(i + 1, n_sites):
                # 跳过已键合的原子对（共价键/配位键不属于碰撞）
                if exclude_bonds and ((i, j) in exclude_bonds or (j, i) in exclude_bonds):
                    continue

                # 使用 pymatgen 的 PBC 感知距离计算
                try:
                    dist = structure.get_distance(i, j)
                except Exception:
                    dist = np.linalg.norm(
                        np.array(structure[i].coords) - np.array(structure[j].coords)
                    )

                vdw_sum = radii[i] + radii[j]
                ratio = dist / vdw_sum if vdw_sum > 0 else float('inf')
                if ratio < min_ratio:
                    min_ratio = ratio

                if dist < self.clash_tolerance * vdw_sum:
                    elem_i = str(structure[i].specie)
                    elem_j = str(structure[j].specie)
                    clashes.append({
                        'i': i,
                        'j': j,
                        'elem_i': elem_i,
                        'elem_j': elem_j,
                        'distance': round(dist, 3),
                        'vdw_sum': round(vdw_sum, 3),
                        'ratio': round(ratio, 3),
                    })

        passed = len(clashes) == 0
        result = {
            'passed': passed,
            'clashes': clashes[:50],  # 最多报告前 50 个碰撞
            'n_clashes': len(clashes),
            'min_ratio': round(min_ratio, 3) if min_ratio != float('inf') else None,
            'tolerance_used': self.clash_tolerance,
        }

        if clashes:
            worst = min(clashes, key=lambda c: c['ratio'])
            result['worst_clash'] = worst
            result['message'] = (
                f"发现 {len(clashes)} 处原子碰撞 "
                f"(最严重: {worst['elem_i']}{worst['i']}—{worst['elem_j']}{worst['j']} "
                f"d={worst['distance']} Å, ratio={worst['ratio']})"
            )
        else:
            result['message'] = f"无碰撞 (最小 vdW ratio={min_ratio:.3f})"

        return result

    # ================================================================
    #  拓扑还原校验
    # ================================================================

    def _check_topology(self, structure, topology_graph, placed_all,
                        bonds_all, lattice_matrix):
        """从 3D 坐标重新推导连通图，与原始拓扑比较。

        方法:
          1. 根据 placed_all 的 bbcode 和 bbtype 将原子分组到 building block
          2. 只考虑 node 类型 building block 之间的连接
             (edge 块作为中间连接器，将连接它的两个 node 块连接起来)
          3. 构建简化图并与 topology_graph 比较
        """
        if placed_all is None or len(placed_all) == 0:
            return {
                'passed': True,
                'skipped': True,
                'message': 'placed_all 为空，跳过拓扑校验'
            }

        # placed_all[i] = [atom_name, x, y, z, charge, cif_label, bbcode, bbtype]
        atom_to_bb = {}
        bbtype_map = {}   # bbcode → {'node', 'edge'}
        for entry in placed_all:
            atom_name = str(entry[0])
            bbcode = int(entry[6])
            bbtype = str(entry[7])
            atom_to_bb[atom_name] = bbcode
            if bbcode not in bbtype_map:
                bbtype_map[bbcode] = bbtype

        # ── 从 bond 数据构建 block 间连接 ──
        # 使用 set 存储邻接关系：{bbcode: set(neighbor_bbcodes)}
        bb_adj = defaultdict(set)
        if bonds_all:
            for bond in bonds_all:
                a1, a2 = str(bond[0]), str(bond[1])
                bb1 = atom_to_bb.get(a1)
                bb2 = atom_to_bb.get(a2)
                if bb1 is not None and bb2 is not None and bb1 != bb2:
                    bb_adj[bb1].add(bb2)
                    bb_adj[bb2].add(bb1)

        # ── 提取 node 类型块之间的连接 ──
        # 对于 node↔edge↔node 的路径，将中间的 edge 替换为直接的 node↔node 边
        node_edges = set()
        node_nodes = set()
        edge_edges = set()

        for bb in bb_adj:
            if bbtype_map.get(bb) == 'node':
                node_nodes.add(bb)
                for nb in bb_adj[bb]:
                    if bbtype_map.get(nb) == 'node':
                        # 直接 node↔node 连接
                        node_edges.add(tuple(sorted((bb, nb))))
                    else:
                        # edge 块：通过其另一端的 node 块间接连接
                        for nb2 in bb_adj[nb]:
                            if nb2 != bb and bbtype_map.get(nb2) == 'node':
                                node_edges.add(tuple(sorted((bb, nb2))))

        # ── 构建简化图 ──
        assembled_graph = nx.Graph()
        assembled_graph.add_nodes_from(node_nodes)
        assembled_graph.add_edges_from(node_edges)

        # ── 与原始拓扑比较 ──
        # TG 是 nx.MultiGraph，转为简单 Graph 比较
        tg_simple = nx.Graph()
        tg_simple.add_nodes_from(topology_graph.nodes())
        tg_simple.add_edges_from(topology_graph.edges())

        n_assembled = assembled_graph.number_of_nodes()
        n_tg = tg_simple.number_of_nodes()
        e_assembled = assembled_graph.number_of_edges()
        e_tg = tg_simple.number_of_edges()

        checks = {
            'n_nodes_match': n_assembled == n_tg,
            'n_edges_match': e_assembled == e_tg,
        }

        # 度序列比较
        assembled_degrees = sorted(d for _, d in assembled_graph.degree())
        tg_degrees = sorted(d for _, d in tg_simple.degree())
        checks['degree_sequence_match'] = (assembled_degrees == tg_degrees)

        # 图同构检验
        isomorphic = False
        if checks['n_nodes_match'] and checks['n_edges_match']:
            try:
                isomorphic = nx.is_isomorphic(assembled_graph, tg_simple)
            except Exception:
                isomorphic = False
        checks['is_isomorphic'] = isomorphic

        all_passed = all(checks.values())

        message_parts = [
            f"节点数: {n_assembled}/{n_tg}",
            f"边数: {e_assembled}/{e_tg}",
            f"同构: {'是' if isomorphic else '否'}",
        ]

        if not checks['n_nodes_match']:
            message_parts.append("节点数不匹配!")
        if not checks['n_edges_match']:
            message_parts.append("边数不匹配!")
        if not checks['degree_sequence_match']:
            message_parts.append("度序列不匹配!")
        if not isomorphic and all_passed is False:
            message_parts.append("图非同构!")

        return {
            'passed': all_passed,
            'n_assembled': n_assembled,
            'n_tg': n_tg,
            'e_assembled': e_assembled,
            'e_tg': e_tg,
            'checks': checks,
            'message': '; '.join(message_parts),
        }

    # ================================================================
    #  对称性守恒校验
    # ================================================================

    def _check_symmetry(self, structure, topology_graph=None):
        """检测生成结构的实际空间群，评估对称性是否严重偏离预期。

        策略:
          - 使用 pymatgen SpacegroupAnalyzer 检测空间群
          - 计算拓扑图的自同构群规模作为"预期对称性"的代理
          - P1 (三斜无对称性) 在拓扑有对称性时视为可疑
        """
        detected_sg = "P1"
        detected_sg_num = 1
        n_symops = 1  # P1 只有恒等操作

        try:
            sga = SpacegroupAnalyzer(structure, symprec=self.symprec)
            detected_sg = sga.get_space_group_symbol()
            detected_sg_num = sga.get_space_group_number()
            # 估算独立对称操作数
            try:
                data = sga.get_symmetry_dataset()
                n_symops = data['n_operations'] if data else 1
            except Exception:
                n_symops = 1
        except Exception as e:
            return {
                'passed': True,
                'skipped': True,
                'detected_sg': 'unknown',
                'message': f"对称性分析异常: {e}",
            }

        # 从拓扑自同构群推断预期对称性水平
        expected_symmetry_high = False
        automorphism_size = 0
        if topology_graph is not None:
            try:
                from symmetry import SymmetryEngine
                # 构建简化的 graph_data 以复用 SymmetryEngine
                n_nodes = len(topology_graph.nodes())
                edges = list(topology_graph.edges())
                simple_graph_data = {
                    'n_nodes': n_nodes,
                    'n_edges': len(edges),
                    'edges': [{'source': int(u), 'target': int(v)}
                              for u, v in edges if str(u).isdigit() and str(v).isdigit()],
                }
                # 仅当 topology_graph 的节点是整数时才可用
                if all(str(n).isdigit() for n in topology_graph.nodes()):
                    engine = SymmetryEngine(simple_graph_data, verbose=False)
                    perms = engine.compute_automorphism_group()
                    automorphism_size = len(perms)
                    expected_symmetry_high = automorphism_size > 1
            except Exception:
                pass

        # 判定标准:
        # 如果拓扑有高对称性但结构是 P1 → 警告
        sg_is_trivial = detected_sg == "P1"
        symmetry_concern = sg_is_trivial and expected_symmetry_high

        # 计算对称性保持度的定性评估
        if sg_is_trivial:
            quality = "poor (P1)" if expected_symmetry_high else "minimal (P1 expected)"
        elif n_symops <= 4:
            quality = "low"
        elif n_symops <= 16:
            quality = "moderate"
        else:
            quality = "high"

        passed = not symmetry_concern

        result = {
            'passed': passed,
            'detected_sg': detected_sg,
            'detected_sg_num': detected_sg_num,
            'n_symops': n_symops,
            'automorphism_size': automorphism_size,
            'symmetry_quality': quality,
        }

        if symmetry_concern:
            result['message'] = (
                f"空间群 {detected_sg} (对称性低), "
                f"但拓扑自同构群有 {automorphism_size} 个元素, "
                f"预期存在非平凡对称性"
            )
        else:
            result['message'] = (
                f"空间群 {detected_sg} (#{detected_sg_num}), "
                f"对称性质量: {quality}"
            )

        return result
