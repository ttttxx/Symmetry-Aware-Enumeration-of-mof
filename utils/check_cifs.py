#mof/utils/check_cifs.py
import os
import sys

from experiment.structure_selection import project_root

# 确保能导入 external_lib
sys.path.append(os.getcwd())
from external_lib.utils.SUPER_check import process_cif


def batch_check(directory):
    if not os.path.exists(directory):
        print(f"Directory not found: {directory}")
        return

    cif_files = [f for f in os.listdir(directory) if f.endswith('.cif')]
    print(f"Checking {len(cif_files)} files in {directory}...")

    passed_count = 0
    for cif in cif_files:
        full_path = os.path.join(directory, cif)
        # process_cif 返回 True 表示通过，False 表示失败
        # 注意：如果失败，它可能会删除文件！建议先备份数据测试。
        try:
            is_valid = process_cif(full_path)
            if is_valid:
                passed_count += 1
            else:
                print(f"  [Failed] {cif} has unreasonable geometry.")
        except Exception as e:
            print(f"  [Error] Could not process {cif}: {e}")

    print(f"Check finished. {passed_count}/{len(cif_files)} passed.")


if __name__ == "__main__":

    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    sys.path.append(project_root)

    # 修改为你刚刚生成的目录路径
    target_dir = os.path.join(project_root,"results/fcu_3c_B_Ch6c_Al_1_Ch_cifs")
    batch_check(target_dir)