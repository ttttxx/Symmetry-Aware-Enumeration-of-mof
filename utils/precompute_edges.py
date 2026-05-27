import os
import json
import sys
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
import warnings

# 忽略 Pymatgen 的繁琐警告
warnings.filterwarnings("ignore")

current_dir = os.path.dirname(os.path.abspath(__file__))
# 获取项目根目录 (mof)
project_root = os.path.dirname(current_dir)
sys.path.append(project_root)


def get_linker_length(cif_path):
    try:
        # 使用 CifParser 读取，确保 primitive=False 以保持原子数一致
        parser = CifParser(cif_path)
        structure = parser.get_structures(primitive=False)[0]

        # 1. 尝试直接从结构中通过 Specie 寻找 (兼容如果用户修改了 CIF 类型为 X 的情况)
        dummy_sites = [site for site in structure if site.specie.symbol == "X"]
        if len(dummy_sites) >= 2:
            p1 = dummy_sites[0].coords
            p2 = dummy_sites[1].coords
            return np.linalg.norm(p1 - p2)

        # 2. 如果方法1失败（即原子类型是 C 但标签是 X），则读取原始标签
        cif_dict = parser.as_dict()
        # 获取第一个数据块的键名
        first_key = list(cif_dict.keys())[0]

        labels = []
        if '_atom_site_label' in cif_dict[first_key]:
            labels = cif_dict[first_key]['_atom_site_label']
            # 确保 labels 是列表（如果只有一个原子，pymatgen有时返回字符串）
            if isinstance(labels, str):
                labels = [labels]

        # 寻找标签以 'X' 开头的原子索引
        dummy_indices = [i for i, label in enumerate(labels) if label.strip().upper().startswith('X')]

        if len(dummy_indices) >= 2:
            # Pymatgen CifParser 通常保持原子顺序与 CIF 文件一致
            # 如果结构中的原子数与标签数匹配，我们可以直接使用索引
            if len(structure) == len(labels):
                p1 = structure.sites[dummy_indices[0]].coords
                p2 = structure.sites[dummy_indices[1]].coords
                return np.linalg.norm(p1 - p2)
            else:
                print(f"Warning: Structure size {len(structure)} != labels size {len(labels)} for {cif_path}")
                return None
        else:
            # print(f"Skipping {cif_path}: Not enough dummy atoms found.")
            return None

    except Exception as e:
        print(f"Error reading {cif_path}: {e}")
        return None


def build_edge_library(edge_dir=os.path.join(project_root, "data/edges"),
                       output_json=os.path.join(project_root, "data/edges/edge_lengths.json")):
    library = {}
    print(f"Scanning {edge_dir}...")

    if not os.path.exists(edge_dir):
        print(f"Error: Directory not found: {edge_dir}")
        return

    files = [f for f in os.listdir(edge_dir) if f.endswith(".cif")]
    print(f"Found {len(files)} CIF files.")

    for f in files:
        full_path = os.path.join(edge_dir, f)
        length = get_linker_length(full_path)
        if length:
            # 名字通常去掉 .cif 后缀
            library[f.replace(".cif", "")] = length

    # 按长度排序
    sorted_lib = dict(sorted(library.items(), key=lambda item: item[1]))

    # 确保输出目录存在
    os.makedirs(os.path.dirname(output_json), exist_ok=True)

    with open(output_json, "w") as f:
        json.dump(sorted_lib, f, indent=4)
    print(f"Saved {len(sorted_lib)} edges to {output_json}")


if __name__ == "__main__":
    build_edge_library()