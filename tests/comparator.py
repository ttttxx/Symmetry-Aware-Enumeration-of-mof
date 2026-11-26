"""
CIF 结构对比模块。

对比参考 CIF 与生成 CIF 的结构差异，返回结构化报告。
"""

import traceback

from pymatgen.core import Structure
from pymatgen.analysis.bond_valence import BVAnalyzer
import numpy as np


class StructureReport:
    """结构化对比报告。"""

    STATUS_OK = "✅正常"
    STATUS_WARN = "⚠️警告"
    STATUS_BAD = "❌严重"

    def __init__(self):
        self.checks = []      # list of (label, status, detail)
        self.warnings = []
        self.errors = []

    def add(self, label, status, detail=""):
        self.checks.append((label, status, detail))
        if status == self.STATUS_WARN:
            self.warnings.append(label)
        elif status == self.STATUS_BAD:
            self.errors.append(label)

    @property
    def has_errors(self):
        return len(self.errors) > 0

    @property
    def has_warnings(self):
        return len(self.warnings) > 0

    def summary(self):
        lines = []
        for label, status, detail in self.checks:
            line = f"  {status} {label}"
            if detail:
                line += f": {detail}"
            lines.append(line)
        return "\n".join(lines)


def compare_structures(ref_path, gen_path):
    """
    对比参考 CIF 与生成 CIF。

    返回 StructureReport 对象。
    """
    report = StructureReport()

    # ── 加载结构 ──
    try:
        ref = Structure.from_file(ref_path)
    except Exception as e:
        report.add("加载参考CIF", StructureReport.STATUS_BAD, traceback.format_exc())
        return report

    try:
        gen = Structure.from_file(gen_path)
    except Exception as e:
        report.add("加载生成CIF", StructureReport.STATUS_BAD, traceback.format_exc())
        return report

    # ── 总原子数 ──
    ref_n = len(ref)
    gen_n = len(gen)
    if ref_n == 0 and gen_n == 0:
        report.add("总原子数", StructureReport.STATUS_BAD, "两者均为空结构")
    elif ref_n == 0:
        report.add("总原子数", StructureReport.STATUS_BAD, f"参考结构原子数为0，生成={gen_n}")
    elif gen_n == 0:
        report.add("总原子数", StructureReport.STATUS_BAD, f"生成结构原子数为0，参考={ref_n}")
    else:
        ratio = gen_n / ref_n
        if 0.5 <= ratio <= 2.0:
            report.add("总原子数", StructureReport.STATUS_OK,
                       f"参考={ref_n}, 生成={gen_n} (比值={ratio:.2f})")
        elif 0.25 <= ratio <= 4.0:
            report.add("总原子数", StructureReport.STATUS_WARN,
                       f"参考={ref_n}, 生成={gen_n} (比值={ratio:.2f}, 偏差较大)")
        else:
            report.add("总原子数", StructureReport.STATUS_BAD,
                       f"参考={ref_n}, 生成={gen_n} (比值={ratio:.2f}, 严重偏差)")

    # ── 晶胞参数 ──
    ref_lat = ref.lattice
    gen_lat = gen.lattice

    ref_vol = ref_lat.volume
    gen_vol = gen_lat.volume

    # 各参数相对偏差
    params = [
        ("a", ref_lat.a, gen_lat.a),
        ("b", ref_lat.b, gen_lat.b),
        ("c", ref_lat.c, gen_lat.c),
        ("α", ref_lat.alpha, gen_lat.alpha),
        ("β", ref_lat.beta, gen_lat.beta),
        ("γ", ref_lat.gamma, gen_lat.gamma),
    ]

    param_diffs = []
    for name, rv, gv in params:
        if abs(rv) < 1e-9:
            continue
        rel = abs(gv - rv) / abs(rv)
        param_diffs.append((name, rv, gv, rel))

    max_param_rel = max(d[3] for d in param_diffs) if param_diffs else 0
    if max_param_rel < 0.05:
        detail = ", ".join(f"{n}: 参考={rv:.3f}, 生成={gv:.3f}" for n, rv, gv, _ in param_diffs)
        report.add("晶胞参数", StructureReport.STATUS_OK, detail)
    elif max_param_rel < 0.20:
        detail = ", ".join(f"{n}: 参考={rv:.3f}, 生成={gv:.3f}" for n, rv, gv, _ in param_diffs)
        report.add("晶胞参数", StructureReport.STATUS_WARN, detail)
    else:
        detail = ", ".join(f"{n}: 参考={rv:.3f}, 生成={gv:.3f}" for n, rv, gv, _ in param_diffs)
        report.add("晶胞参数", StructureReport.STATUS_BAD, detail)

    # ── 晶胞体积 ──
    if min(ref_vol, gen_vol) > 0:
        vol_ratio = gen_vol / ref_vol
        if vol_ratio > 2.0 or vol_ratio < 0.5:
            report.add("晶胞体积", StructureReport.STATUS_BAD,
                       f"参考={ref_vol:.1f} Å³, 生成={gen_vol:.1f} Å³ (比值={vol_ratio:.2f}, >50%差异)")
        else:
            report.add("晶胞体积", StructureReport.STATUS_OK,
                       f"参考={ref_vol:.1f} Å³, 生成={gen_vol:.1f} Å³ (比值={vol_ratio:.2f})")
    else:
        report.add("晶胞体积", StructureReport.STATUS_BAD, "晶胞体积为零或负值")

    # ── 空间群 ──
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        ref_sg = SpacegroupAnalyzer(ref, symprec=0.3).get_space_group_symbol()
    except Exception:
        ref_sg = "?"
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        gen_sg = SpacegroupAnalyzer(gen, symprec=0.3).get_space_group_symbol()
    except Exception:
        gen_sg = "?"

    if ref_sg and gen_sg:
        if ref_sg == gen_sg:
            report.add("空间群", StructureReport.STATUS_OK, f"{ref_sg}")
        else:
            report.add("空间群", StructureReport.STATUS_WARN,
                       f"参考={ref_sg}, 生成={gen_sg}")
    else:
        report.add("空间群", StructureReport.STATUS_OK, f"参考={ref_sg or 'P1'}, 生成={gen_sg or 'P1'} (默认为P1)")

    # ── 元素组成 ──
    ref_comp = ref.composition.get_el_amt_dict()
    gen_comp = gen.composition.get_el_amt_dict()
    all_elems = set(ref_comp.keys()) | set(gen_comp.keys())

    comp_diffs = []
    for el in sorted(all_elems, key=str):
        rc = ref_comp.get(el, 0)
        gc = gen_comp.get(el, 0)
        if abs(rc) < 1e-9 and abs(gc) < 1e-9:
            continue
        rel = abs(gc - rc) / max(abs(rc), 1)
        comp_diffs.append((str(el), rc, gc, rel))

    max_comp_rel = max(d[3] for d in comp_diffs) if comp_diffs else 0
    if max_comp_rel < 0.10:
        detail = ", ".join(f"{el}: 参考={rc:.1f}, 生成={gc:.1f}" for el, rc, gc, _ in comp_diffs)
        report.add("元素组成", StructureReport.STATUS_OK, detail)
    elif max_comp_rel < 0.30:
        detail = ", ".join(f"{el}: 参考={rc:.1f}, 生成={gc:.1f}" for el, rc, gc, _ in comp_diffs)
        report.add("元素组成", StructureReport.STATUS_WARN, detail)
    else:
        detail = ", ".join(f"{el}: 参考={rc:.1f}, 生成={gc:.1f}" for el, rc, gc, _ in comp_diffs)
        report.add("元素组成", StructureReport.STATUS_BAD, detail)

    # ── 原子间距检查 ──
    ref_min = _min_interatomic_dist(ref)
    gen_min = _min_interatomic_dist(gen)

    # 生成结构的最近邻成键范围（用于判断短间距是否为键内原子）
    gen_nn_max = _max_nearest_neighbor_dist(gen) if gen_n > 1 else None

    if ref_min is not None and ref_min < 0.8:
        report.add("参考结构原子间距", StructureReport.STATUS_BAD,
                   f"最短原子间距={ref_min:.3f}Å (<0.8Å, 原子重叠)")
    elif ref_min is not None:
        report.add("参考结构原子间距", StructureReport.STATUS_OK,
                   f"最短原子间距={ref_min:.3f}Å")

    if gen_min is not None and gen_min < 0.8:
        if gen_nn_max is not None and gen_nn_max > 0.8:
            report.add("生成结构原子间距", StructureReport.STATUS_WARN,
                       f"最短原子间距={gen_min:.3f}Å (<0.8Å)，"
                       f"但最近邻范围包含正常成键值({gen_min:.2f}~{gen_nn_max:.2f}Å)，"
                       f"可能为键内原子，建议人工确认")
        else:
            report.add("生成结构原子间距", StructureReport.STATUS_BAD,
                       f"最短原子间距={gen_min:.3f}Å (<0.8Å, 原子重叠)")
    elif gen_min is not None:
        report.add("生成结构原子间距", StructureReport.STATUS_OK,
                   f"最短原子间距={gen_min:.3f}Å")

    # Bond distance sanity: 检查是否有 >3.0Å 的"键"（可能缺键）
    if gen_n > 1:
        if gen_nn_max is not None and gen_nn_max > 3.0:
            report.add("生成结构成键间距", StructureReport.STATUS_WARN,
                       f"最近邻最大距离={gen_nn_max:.2f}Å (>3Å, 可能存在缺键)")
        elif gen_nn_max is not None:
            report.add("生成结构成键间距", StructureReport.STATUS_OK,
                       f"最近邻距离范围={gen_min:.2f}~{gen_nn_max:.2f}Å")

    return report


def _min_interatomic_dist(structure):
    """计算结构中最短原子间距（排除自身）。"""
    min_d = float("inf")
    n = len(structure)
    if n < 2:
        return min_d
    # 使用距离矩阵（采样前500个原子以避免超大结构O(N²)爆炸）
    indices = list(range(min(n, 500)))
    for i_idx, i in enumerate(indices):
        for j in indices[i_idx + 1:]:
            d = structure.get_distance(i, j)
            if d < min_d:
                min_d = d
    return min_d


def _max_nearest_neighbor_dist(structure):
    """计算结构中每个原子的最近邻距离，返回最大值。"""
    n = len(structure)
    if n < 2:
        return None
    # 使用 pymatgen 的 get_all_neighbors，半径取 5Å 覆盖成键范围
    from pymatgen.analysis.local_env import JmolNN
    try:
        nn_alg = JmolNN()
        max_nn = 0.0
        for i in range(min(n, 500)):
            neighbors = nn_alg.get_nn_info(structure, i)
            if neighbors:
                for nb in neighbors:
                    d = structure.get_distance(i, nb["site_index"])
                    if d > max_nn:
                        max_nn = d
        return max_nn if max_nn > 0 else None
    except Exception:
        # 回退：暴力计算最近邻
        max_min_dist = 0.0
        indices = list(range(min(n, 200)))
        for i_idx, i in enumerate(indices):
            min_dist = float("inf")
            for j in range(n):
                if i != j:
                    d = structure.get_distance(i, j)
                    if d < min_dist:
                        min_dist = d
            if min_dist > max_min_dist and min_dist < 5.0:
                max_min_dist = min_dist
        return max_min_dist if max_min_dist > 0 else None


# ── 便捷函数 ──

def compare_and_print(ref_path, gen_path):
    """对比并打印报告，返回 report 对象。"""
    report = compare_structures(ref_path, gen_path)
    print(report.summary())
    return report
