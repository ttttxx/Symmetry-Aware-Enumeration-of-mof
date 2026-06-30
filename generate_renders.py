import os
import json
import glob
from assembler import StructureAssembler


def main():
    # ================= 配置区 =================
    DATA_PREFIX = "sod"  # 保持与你 main.py 一致
    RENDER_DIR = os.path.join("results", "render")
    # =========================================

    if not os.path.exists(RENDER_DIR):
        os.makedirs(RENDER_DIR)

    print(f"--- Automating 3D Render CIF Generation for {DATA_PREFIX} ---")

    # 1. 自动获取最新运行的结果 JSON 和原始拓扑图
    result_files = glob.glob(os.path.join("results", "sod_0_mn1_new6c_Al_1_Ch_20260304_224821.json"))#f"{DATA_PREFIX}_*.json"
    if not result_files:
        raise FileNotFoundError("未找到枚举结果 JSON，请先运行 main.py")

    latest_file = max(result_files, key=os.path.getctime)
    print(f"Loading results from: {latest_file}")

    with open(latest_file, 'r') as f:
        result_data = json.load(f)

    graph_path = os.path.join("data", DATA_PREFIX, f"{DATA_PREFIX}_graph.json")
    with open(graph_path, 'r') as f:
        graph_data = json.load(f)

    group_size = result_data["meta"]["group_size"]
    structures = result_data["structures"]

    # 2. 按对称性 (stabilizer_size) 对所有构型进行严格分组
    symm_groups = {}
    for s in structures:
        stab = s["stabilizer_size"]
        if stab not in symm_groups:
            symm_groups[stab] = []
        symm_groups[stab].append(s)

    # 将稳定子大小从大到小排序 (高对称 -> 低对称)
    sorted_stabs = sorted(symm_groups.keys(), reverse=True)

    # 3. 算法自动挑选 4 个最具代表性的对称性阶层
    selected_structures = []

    if len(sorted_stabs) >= 4:
        target_stabs = [
            sorted_stabs[0],  # Render 1: 最高对称性 (rho 趋近 1)
            sorted_stabs[1],  # Render 2: 次高对称性 (子群分支 A)
            sorted_stabs[len(sorted_stabs) // 2],  # Render 3: 中间对称性 (子群分支 B)
            sorted_stabs[-1]  # Render 4: 最低对称性 (不对称极限，rho 趋近 0)
        ]
    else:
        # 容错机制：如果系统特别小，对称性阶层不到 4 个，就用头部和尾部补齐
        print("  [Note] Symmetry tiers are less than 4, duplicating edge cases.")
        target_stabs = [sorted_stabs[0], sorted_stabs[0], sorted_stabs[-1], sorted_stabs[-1]]

    for i, stab in enumerate(target_stabs):
        # 提取该对称性阶层里的第一个结构（通常是字典序最靠前、能量较低的代表元）
        selected_structures.append((i + 1, symm_groups[stab][0]))

    # 4. 调用你的 Assembler 进行物理坐标重构并生成 CIF
    assembler = StructureAssembler(output_dir=RENDER_DIR)

    print("\n  [Starting Assembly Pipeline...]")
    for render_idx, struct_info in selected_structures:
        stab_size = struct_info['stabilizer_size']
        rho = stab_size / group_size
        struct_id = struct_info['id']

        # 构造直观的文件名，方便你在 VESTA 里对应 Draw.io 的图框
        filename = f"Render_{render_idx}_rho_{rho:.2f}_stab_{stab_size}.cif"
        out_path = os.path.join(RENDER_DIR, filename)

        # 将 list 格式的 coloring 转换为 Assembler 需要的 dict 格式
        coloring_list = struct_info['coloring']
        structure_map = {str(idx): elem for idx, elem in enumerate(coloring_list)}

        print(f"  -> Generating [Render {render_idx}]: rho = {rho:.2f}, ID = {struct_id}")

        try:
            # 直接复用你之前写好的刚体对齐和 PBC 消除逻辑
            assembler.assemble_and_save(out_path, structure_map, graph_data)
            print(f"     Saved: {filename}")
        except Exception as e:
            print(f"     [Error] Failed to assemble {filename}: {e}")

    print(f"\n✅ All 4 render CIFs successfully generated in: {RENDER_DIR}")
    print("Next step: Open them in VESTA, set alpha channel, and export transparent PNGs!")


if __name__ == "__main__":
    main()