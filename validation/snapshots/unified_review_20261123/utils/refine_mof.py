from ase.io import read, write
from ase.optimize import BFGS
# from xtb.ase.calculator import XTB
import os


def refine_mof_structure(input_cif, output_cif):
    """
    使用 xTB 自动弛豫 MOF 结构并恢复对称性
    """
    print(f"--- Refining {input_cif} ---")

    # 1. 读取结构
    atoms = read(input_cif)

    # 2. 设置计算器 (GFN-FF 是专门为大型分子和 MOF 设计的力场水平，极快)
    # 如果需要更高精度，可以使用 method='GFN2-xTB'
    # atoms.calc = XTB(method='GFN-FF')

    # 3. 运行几何优化 (弛豫原子位置)
    # 注意：对于 MOF，通常需要固定晶格只优化原子位置，或者使用 ExpCellFilter 同时优化晶格
    from ase.constraints import ExpCellFilter
    ecf = ExpCellFilter(atoms)

    dyn = BFGS(ecf, trajectory='relax.traj', logfile='relax.log')
    dyn.run(fmax=0.05)  # 收敛阈值

    # 4. 保存优化后的结构
    write(output_cif, atoms)

    # 5. 再次使用 Pymatgen 恢复对称性 (因为优化后坐标会微动)
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    refined_pymatgen_struct = Structure.from_file(output_cif)
    sga = SpacegroupAnalyzer(refined_pymatgen_struct, symprec=0.1)

    print(f"Final Spacegroup: {sga.get_space_group_symbol()}")
    final_struct = sga.get_refined_structure()
    final_struct.to(filename=output_cif)

# 调用示例
# refine_mof_structure("results/cifs/my_mof.cif", "results/cifs/my_mof_optimized.cif")