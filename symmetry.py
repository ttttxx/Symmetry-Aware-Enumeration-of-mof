# mof/symmetry.py
import numpy as np
import os
import json
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class SymmetryEngine:
    """
    SymmetryEngine (最终修复版):
    1. 修复 get_all_distances 返回二维数组导致的 '全是 Node 0' 映射 bug。
    2. 包含之前所有的增强逻辑 (元素清洗、Lattice 转置、容差放宽)。
    """

    def __init__(self, graph_json_data, verbose=True):
        self.verbose = verbose
        # 确保按 ID 排序
        self.nodes = sorted(graph_json_data["nodes"], key=lambda x: x['id'])
        self.num_nodes = len(self.nodes)

        # --- 修复 1: Lattice 转置 ---
        raw_lattice = np.array(graph_json_data["lattice"])
        lattice_mat = raw_lattice.T

        # --- 修复 2: 元素清洗 ---
        raw_species = [n["element"] for n in self.nodes]
        clean_species = []
        for sp in raw_species:
            if "Z" in sp or "M" in sp:
                clean_species.append("Zn")
            else:
                clean_species.append("C")

        coords = [n.get("frac", n.get("coord")) for n in self.nodes]
        self.structure = Structure(lattice_mat, clean_species, coords)

        if self.verbose:
            print(f"  [Init] Structure loaded: {self.structure.formula} ({self.num_nodes} atoms)")

    def _get_robust_analyzer(self):
        precisions = [0.01, 0.1, 0.3]
        best_sga = None
        best_num_ops = 0

        for prec in precisions:
            try:
                sga = SpacegroupAnalyzer(self.structure, symprec=prec)
                ops = sga.get_symmetry_operations()
                if len(ops) > best_num_ops:
                    best_num_ops = len(ops)
                    best_sga = sga
                if len(ops) > 16: break
            except Exception:
                continue
        return best_sga

    def compute_automorphism_group(self):
        if self.verbose:
            print(f"  [Symmetry] Analyzing symmetry...")

        sga = self._get_robust_analyzer()
        if sga is None:
            return [tuple(range(self.num_nodes))]

        symm_ops = sga.get_symmetry_operations()
        permutations = []
        frac_coords = self.structure.frac_coords

        # 容差
        DIST_TOL = 0.8

        for op in symm_ops:
            # 过滤镜像 (保留纯旋转)
            if np.linalg.det(op.rotation_matrix) < 0:
                continue

            current_perm = []
            is_valid_perm = True

            for i in range(self.num_nodes):
                original_coord = frac_coords[i]
                new_coord = op.operate(original_coord)

                # --- 核心修复：添加 .flatten() ---
                # Pymatgen 可能返回 [[d1, d2...]] (shape 1xN)，这会导致 np.where 返回行索引 0
                # flatten() 确保它变回 [d1, d2...] (shape N)
                dists = self.structure.lattice.get_all_distances(new_coord, frac_coords).flatten()

                # 寻找匹配
                match_indices = np.where(dists < DIST_TOL)[0]

                if len(match_indices) == 0:
                    is_valid_perm = False
                    break

                # 只要是一维数组，这里取出的就是正确的节点 ID
                target_idx = match_indices[0]
                current_perm.append(target_idx)

            if is_valid_perm:
                if len(set(current_perm)) == self.num_nodes:
                    permutations.append(tuple(current_perm))

        # 强制包含恒等变换
        identity = tuple(range(self.num_nodes))
        if identity not in permutations:
            permutations.append(identity)

        unique_perms = sorted(list(set(permutations)))

        if self.verbose:
            print(f"  [Symmetry] Generated {len(unique_perms)} valid permutations (from {len(symm_ops)} ops).")

        return unique_perms


if __name__ == "__main__":
    # 简单的自测代码
    import sys


    json_path = "data/sod/sod_graph.json"


    if json_path:
        print(f"--- Testing fix on {json_path} ---")
        with open(json_path, 'r') as f:
            data = json.load(f)
        engine = SymmetryEngine(data, verbose=True)
        perms = engine.compute_automorphism_group()
        print(f"Final Group Size: {len(perms)}")
        if len(perms) > 1:
            print("[SUCCESS] Found symmetry!")
            print(f"Sample perm: {perms[1]}")
        else:
            print("[WARNING] Still only found Identity.")
    else:
        print("Test file not found.")