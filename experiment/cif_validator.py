import os
import glob
from pymatgen.core import Structure
from chgnet.model import CHGNet
import pandas as pd


def validate_cifs(cif_directory):
    """
    读取指定目录下的所有 CIF 文件，使用 CHGNet 预测其热力学稳定性（能量）。
    """
    print("Loading Universal ML Potential (CHGNet)...")
    model = CHGNet.load()  # 首次运行会自动下载预训练权重

    cif_files = glob.glob(os.path.join(cif_directory, "*.cif"))
    results = []

    print(f"Found {len(cif_files)} CIF files. Starting thermodynamic evaluation...")

    for cif_path in cif_files:
        filename = os.path.basename(cif_path)
        try:
            # 1. 读取 CIF 结构
            struct = Structure.from_file(cif_path)

            # 2. 预测能量 (无需进行耗时的全弛豫，直接计算单点能即可快速筛选)
            prediction = model.predict_structure(struct)

            # 提取预测的每原子能量 (eV/atom)
            energy_per_atom = prediction['e']

            # 可选：如果你也想验证力学稳定性（受力是否趋近于0）
            max_force = max([abs(f) for force in prediction['f'] for f in force])

            results.append({
                "Filename": filename,
                "Energy_per_atom_eV": energy_per_atom,
                "Max_Force_eV_A": max_force,
                "Status": "Success"
            })
            print(f"Evaluated {filename}: Energy = {energy_per_atom:.4f} eV/atom")

        except Exception as e:
            print(f"Error evaluating {filename}: {e}")
            results.append({
                "Filename": filename,
                "Energy_per_atom_eV": None,
                "Max_Force_eV_A": None,
                "Status": f"Failed: {str(e)}"
            })

    # 将结果保存为 CSV，用于替换你论文图 4 (Fig 4) 的数据源
    df = pd.DataFrame(results)
    output_csv = "thermodynamic_validation_results.csv"
    df.to_csv(output_csv, index=False)
    print(f"\nValidation complete. Real physical energies saved to {output_csv}")


if __name__ == "__main__":
    # 指定你的生成 CIF 文件所在的目录，例如：
    CIF_DIR = "data/mof"
    validate_cifs(CIF_DIR)