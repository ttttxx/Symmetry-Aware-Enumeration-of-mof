import os
from pymatgen.core import Structure
# TODO：需不需要加入main主流程

def check_overlaps(directory, cutoff=0.8):
    """
    检查指定目录下所有 CIF 是否存在原子间距小于 cutoff 的情况。
    cutoff: 阈值，单位埃。一般金属-有机键不应小于 1.0-1.5A，此处 0.8 是极保守的防重叠阈值。
    """
    cif_files = [f for f in os.listdir(directory) if f.endswith('.cif')]
    print(f"Checking {len(cif_files)} files in {directory} for overlaps < {cutoff} A...")

    for cif in cif_files:
        path = os.path.join(directory, cif)
        try:
            struct = Structure.from_file(path)
            # 获取结构中所有原子对的距离矩阵
            dist_matrix = struct.distance_matrix
            # 将对角线（自己到自己）设为无穷大，避免误判
            num_sites = len(struct)
            for i in range(num_sites):
                dist_matrix[i, i] = float('inf')

            min_dist = dist_matrix.min()

            if min_dist < cutoff:
                print(f"  [WARNING] {cif}: Clash detected! Min distance = {min_dist:.3f} A")
            else:
                # 还可以检查密度是否异常，作为辅助判断
                # print(f"  [Pass] {cif} (Min dist: {min_dist:.3f} A)")
                pass

        except Exception as e:
            print(f"  [Error] Cannot read {cif}: {e}")


if __name__ == "__main__":
    target_dir = "results/sod_6c_CdZn_1_Ch8c_Co_1_Ch_cifs"
    check_overlaps(target_dir)