#!/usr/bin/env python
"""
MOF 结构组装测试主入口。

流程：
  1. 扫描 tests/tobacco_output_cifs/ 下所有参考 CIF
  2. 交互选择测试范围
  3. 解析文件名 → 调用组装管线生成 CIF
  4. 对比参考与生成 CIF，输出汇总报告
"""

import os
import sys
import glob
import random
import shutil
from collections import defaultdict
from datetime import datetime

# 确保项目根目录在 sys.path 中
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT_ROOT)
# 【新增此行】将当前工作目录强制切换到项目根目录，
# 确保底层通过相对路径 "data/..." 查找拓扑文件的模块不会报错
os.chdir(PROJECT_ROOT)

from tests.name_parser import parse_filename
from tests.comparator import compare_structures, StructureReport


def scan_reference_cifs():
    """扫描参考 CIF 目录，返回排序后的文件名列表。"""
    ref_dir = os.path.join(os.path.dirname(__file__), "tobacco_output_cifs")
    if not os.path.isdir(ref_dir):
        print(f"[错误] 参考 CIF 目录不存在: {ref_dir}")
        return []
    files = sorted(
        f for f in os.listdir(ref_dir) if f.endswith(".cif")
    )
    return files


def select_files(all_files):
    """交互选择待测试的文件。返回选中文件名列表。"""
    print(f"\n发现 {len(all_files)} 个参考 CIF 文件：")
    for i, f in enumerate(all_files, 1):
        print(f"  [{i}] {f}")

    print("\n选择模式：")
    print("  [1] 手动输入编号（逗号分隔，如 1,3,5）")
    print("  [2] 随机抽样 N 个")
    print("  [3] 全部")

    while True:
        choice = input("\n请输入选择 [1/2/3]: ").strip()
        if choice == "1":
            while True:
                raw = input("请输入编号（逗号分隔）: ").strip()
                try:
                    indices = [int(x.strip()) - 1 for x in raw.split(",")]
                    selected = []
                    for idx in indices:
                        if 0 <= idx < len(all_files):
                            selected.append(all_files[idx])
                        else:
                            print(f"  编号 {idx + 1} 超出范围，已忽略")
                    if selected:
                        return selected
                    print("  未选中有效文件，请重新输入")
                except ValueError:
                    print("  格式错误，请输入如 1,3,5")
        elif choice == "2":
            while True:
                raw = input("随机抽取几个？ ").strip()
                try:
                    n = int(raw)
                    if 1 <= n <= len(all_files):
                        return random.sample(all_files, n)
                    print(f"  请输入 1~{len(all_files)} 之间的数字")
                except ValueError:
                    print("  请输入有效数字")
        elif choice == "3":
            return list(all_files)
        else:
            print("  请输入 1、2 或 3")


def run_assembly_for_case(ref_filename):
    """
    对单个参考 CIF 文件名运行组装管线。

    返回 (generated_cif_path, error_msg)。
    如果组装成功，error_msg 为 None；否则 generated_cif_path 为 None。
    """
    tmp_dir = os.path.join(os.path.dirname(__file__), "tmp_output")
    os.makedirs(tmp_dir, exist_ok=True)

    # 1. 解析文件名
    try:
        topology, nodes, edge = parse_filename(ref_filename)
    except Exception as e:
        return None, f"文件名解析失败: {e}"

    print(f"  拓扑={topology}, 节点={nodes}, edge={edge}")

    # 2. 检查拓扑数据是否存在
    graph_path = os.path.join(PROJECT_ROOT, "data", topology, f"{topology}_graph.json")
    if not os.path.exists(graph_path):
        return None, f"拓扑数据不存在: data/{topology}/ (跳过)"

    # 3. 构建 user_selected_ligands
    user_selected_ligands = defaultdict(list)
    for node_name in nodes:
        # 从 "4c_Cu_1_Ch" 中提取 connectivity degree
        try:
            deg = int(node_name.split("c")[0])
        except (ValueError, IndexError):
            return None, f"无法从节点名 '{node_name}' 解析连通度"
        if node_name not in user_selected_ligands[deg]:
            user_selected_ligands[deg].append(node_name)
    user_selected_ligands = dict(user_selected_ligands)

    # 4. 调用组装管线
    try:
        from main import run_enumeration
        cif_paths, _, _ = run_enumeration(
            data_prefix=topology,
            user_selected_ligands=user_selected_ligands,
            forced_edge=edge,
            max_cifs=1,
            max_results=5000,
        )
    except Exception as e:
        return None, f"组装管线异常: {e}"

    if not cif_paths:
        return None, "组装管线未生成任何 CIF（可能因 shape mismatch 或验证失败）"

    # 5. 复制到 tmp_output
    src = cif_paths[0]
    base = os.path.splitext(ref_filename)[0]
    dst = os.path.join(tmp_dir, f"{base}_gen.cif")
    shutil.copy2(src, dst)
    print(f"  生成 CIF → {dst}")

    return dst, None


def merge_reports(reports):
    """
    合并多个报告的相同错误。
    返回 (merged_lines, bug_groups)。

    bug_groups: {bug_signature: [case_names]}
    """
    bug_groups = defaultdict(list)

    for case_name, report in reports:
        # 用 (error_labels_tuple) 作为 bug 签名
        if report.errors:
            sig = tuple(sorted(report.errors))
        elif report.warnings:
            sig = tuple(sorted(report.warnings))
        else:
            sig = ("passed",)
        bug_groups[sig].append(case_name)

    lines = []
    lines.append("=" * 60)
    lines.append("汇总报告")
    lines.append("=" * 60)

    total = len(reports)
    passed = sum(1 for _, r in reports if not r.has_errors)
    lines.append(f"\n总计: {total} 个用例, 通过: {passed}, 失败: {total - passed}")

    # 按 bug 分组输出
    bug_idx = 1
    for sig, cases in sorted(bug_groups.items(), key=lambda x: -len(x[1])):
        if sig == ("passed",):
            continue
        lines.append(f"\n{'─' * 60}")
        lines.append(f"Bug #{bug_idx} — 影响 {len(cases)} 个用例:")
        bug_idx += 1
        for c in cases:
            lines.append(f"  • {c}")
        # 取第一个用例的详细报告
        first_report = next(r for cn, r in reports if cn == cases[0])
        lines.append(f"\n  代表用例 [{cases[0]}] 详细报告:")
        lines.append(first_report.summary())

    # 通过的用例
    passed_cases = [cn for cn, r in reports if not r.has_errors]
    if passed_cases:
        lines.append(f"\n{'─' * 60}")
        lines.append(f"通过 ({len(passed_cases)} 个):")
        for c in passed_cases:
            lines.append(f"  ✅ {c}")

    return "\n".join(lines), dict(bug_groups)


# ── 问题分析与修复阶段 ──

def _get_function_signatures(module_path):
    """读取模块文件，只提取函数/类/方法签名（不读完整实现）。"""
    sigs = []
    if not os.path.exists(module_path):
        return sigs
    with open(module_path, "r") as f:
        for line in f:
            stripped = line.rstrip()
            if stripped.startswith("class ") or stripped.startswith("    def ") or stripped.startswith("def "):
                sigs.append(stripped)
    return sigs


def _infer_assembly_cause(error_msg):
    """根据组装错误信息推断可能原因。"""
    if "拓扑数据不存在" in error_msg:
        return {
            "category": "路径/配置问题",
            "cause": "参考 CIF 的拓扑名在 data/ 目录下没有对应的 graph JSON 文件。",
            "suggestion": "检查是否缺少拓扑数据文件，或拓扑名解析有误。",
        }
    if "文件名解析失败" in error_msg:
        return {
            "category": "文件名格式问题",
            "cause": "CIF 文件名不符合预期的 {topo}_v1-{node}_..._{edge}.cif 格式。",
            "suggestion": "检查文件名拼写，确认 vN- 节点标记和 edge 名是否正确。",
        }
    if "连通度" in error_msg:
        return {
            "category": "参数问题",
            "cause": "节点名中无法提取连通度（如 4c_xxx 中的数字+c 模式）。",
            "suggestion": "确认节点 CIF 文件名格式为 {N}c_{name}.cif。",
        }
    if "组装管线异常" in error_msg:
        return {
            "category": "代码 bug",
            "cause": f"组装过程中抛出未捕获异常: {error_msg}",
            "suggestion": "需检查对应模块的堆栈信息，定位具体出错函数。",
        }
    if "未生成任何 CIF" in error_msg:
        return {
            "category": "验证/匹配失败",
            "cause": "枚举完成但组装未产出 CIF（shape mismatch / 验证失败 / 无合法着色方案）。",
            "suggestion": "检查 topology_requirements 与所选配体是否匹配，或放宽验证条件。",
        }
    return {
        "category": "未知",
        "cause": error_msg,
        "suggestion": "请人工检查错误信息。",
    }


def _map_check_to_module(check_label):
    """将对比报告中的检查标签映射到可能的代码模块和函数。"""
    mapping = {
        "总原子数": {
            "module": "assembler.py",
            "functions": ["StructureAssembler.assemble_and_save", "calculate_X_X_bonds_and_update_charges"],
            "note": "原子放置/融合/去重逻辑可能导致原子数偏差",
        },
        "晶胞参数": {
            "module": "assembler.py / main.py",
            "functions": ["StructureAssembler._compute_lattice_scale", "enrich_graph_data_with_topology"],
            "note": "晶胞缩放或拓扑矩阵计算可能引入晶胞参数偏差",
        },
        "晶胞体积": {
            "module": "assembler.py",
            "functions": ["StructureAssembler._compute_lattice_scale", "StructureAssembler.assemble_and_save"],
            "note": "晶胞缩放因子计算或 unit cell 构造可能有问题",
        },
        "空间群": {
            "module": "assembler.py",
            "functions": ["StructureAssembler.assemble_and_save (SpacegroupAnalyzer 部分)"],
            "note": "结构对称化步骤或原子位置精度可能影响空间群识别",
        },
        "元素组成": {
            "module": "assembler.py",
            "functions": ["StructureAssembler.assemble_and_save (merge/fusion 部分)"],
            "note": "X 原子融合、merge_sites 或原子去重可能改变元素计数",
        },
        "参考结构原子间距": {
            "module": "N/A (参考 CIF 自身问题)",
            "functions": [],
            "note": "参考 CIF 本身存在原子重叠，非本组装代码问题",
        },
        "生成结构原子间距": {
            "module": "assembler.py",
            "functions": ["get_template_topology_info", "adjust_placed_edges_to_fit_periodicity"],
            "note": "原子放置位置可能导致重叠或过近",
        },
        "生成结构成键间距": {
            "module": "assembler.py",
            "functions": ["get_template_topology_info", "calculate_X_X_bonds_and_update_charges"],
            "note": "键连接构建或 X-X 键计算可能遗漏成键",
        },
    }
    return mapping.get(check_label, {
        "module": "待定位",
        "functions": [],
        "note": "需要人工分析",
    })


def _collect_problems(case_details, bug_groups):
    """
    从 case_details 和 bug_groups 中收集所有问题，去重后返回扁平列表。
    每项 dict 包含：idx(填充时), case_name, ref_file, affected_count,
                   category, error, report, label, status, detail, diag
    """
    problems = []
    seen_sigs = set()

    for sig, case_names in sorted(bug_groups.items(), key=lambda x: -len(x[1])):
        if sig == ("passed",):
            continue

        sig_key = tuple(sorted(sig))
        if sig_key in seen_sigs:
            continue
        seen_sigs.add(sig_key)

        rep_case = case_names[0]
        rep_detail = next(cd for cd in case_details if cd["case_name"] == rep_case)
        error = rep_detail.get("error")
        report = rep_detail.get("report")

        if error:
            diag = _infer_assembly_cause(error)
            problems.append({
                "case_name": rep_case,
                "ref_file": rep_detail["ref_file"],
                "affected_count": len(case_names),
                "category": "assembly",
                "diag": diag,
                "error": error,
                "report": report,
                "label": f"组装报错 — {diag['category']}",
                "status": StructureReport.STATUS_WARN,
                "detail": error,
                "module": "",
                "functions": [],
                "note": diag.get("cause", ""),
            })
        elif report:
            failing_checks = [
                (label, status, detail)
                for label, status, detail in report.checks
                if status in (StructureReport.STATUS_BAD, StructureReport.STATUS_WARN)
            ]
            for label, status, detail in failing_checks:
                mod_info = _map_check_to_module(label)
                problems.append({
                    "case_name": rep_case,
                    "ref_file": rep_detail["ref_file"],
                    "affected_count": len(case_names),
                    "category": "comparison",
                    "report": report,
                    "error": None,
                    "label": label,
                    "status": status,
                    "detail": detail,
                    "module": mod_info["module"],
                    "functions": mod_info["functions"],
                    "note": mod_info["note"],
                })

    return problems


def _select_problems(problems):
    """列出问题并让用户选择需要修复的编号。返回选中子列表。"""
    print(f"\n{'─' * 60}")
    print(f"发现 {len(problems)} 类问题：\n")
    for i, p in enumerate(problems, 1):
        print(f"  [{i}] {p['status']} {p['label']}")
        print(f"      用例: {p['case_name']} (影响 {p['affected_count']} 个)")
        print(f"      详情: {p['detail']}")

    while True:
        raw = input("\n请选择需要修复的问题编号（逗号分隔，回车跳过）: ").strip()
        if not raw:
            return []
        try:
            indices = [int(x.strip()) - 1 for x in raw.split(",")]
            selected = []
            for idx in indices:
                if 0 <= idx < len(problems):
                    selected.append(problems[idx])
                else:
                    print(f"  编号 {idx + 1} 超出范围，已忽略")
            if selected:
                return selected
            print("  未选中有效问题，请重新输入")
        except ValueError:
            print("  格式错误，请输入如 1,3,5")


def _retest_single_case(problem):
    """重测单个用例：重新组装 + 对比，返回新的 StructureReport 或 None。"""
    ref_file = problem.get("ref_file", "")
    if not ref_file:
        print("    [重测] 缺少 ref_file，无法重测")
        return None

    ref_path = os.path.join(os.path.dirname(__file__), "tobacco_output_cifs", ref_file)
    gen_path, error = run_assembly_for_case(ref_file)

    if error:
        print(f"    [重测] 组装失败: {error}")
        return None

    return compare_structures(ref_path, gen_path)


def _interactive_fix_loop(case_details, bug_groups):
    """交互式修复循环：收集问题 → 用户选择 → 调用 RepairAgent 逐个修复。"""
    print("\n" + "=" * 60)
    print("问题分析与修复阶段")
    print("=" * 60)

    has_failures = any(
        sig != ("passed",)
        for sig in bug_groups
    )
    if not has_failures:
        print("\n所有用例均通过，无需修复。")
        return

    total_failed = sum(
        len(cases) for sig, cases in bug_groups.items() if sig != ("passed",)
    )
    print(f"\n共 {len(bug_groups) - (1 if ('passed',) in bug_groups else 0)} 类错误，影响 {total_failed} 个用例。")

    while True:
        ans = input("\n是否进入自动修复流程？[y/n]: ").strip().lower()
        if ans == "y":
            break
        elif ans == "n":
            print("跳过修复流程。")
            return
        print("  请输入 y 或 n")

    # 1. 收集所有问题（同类去重）
    problems = _collect_problems(case_details, bug_groups)

    # 2. 用户选择
    selected = _select_problems(problems)
    if not selected:
        print("未选择任何问题，退出修复流程。")
        return

    # 3. 初始化 RepairAgent
    try:
        from tests.repair_agent import RepairAgent
        agent = RepairAgent()
    except RuntimeError as e:
        print(f"\n无法初始化 RepairAgent: {e}")
        print("请设置环境变量 OPENAI_API_KEY 后重试。")
        return
    except Exception as e:
        print(f"\n导入 RepairAgent 失败: {e}")
        return

    # 4. 逐个修复
    fixed_count = 0
    for i, problem in enumerate(selected, 1):
        print(f"\n{'#' * 60}")
        print(f"修复进度: [{i}/{len(selected)}]")
        print(f"{'#' * 60}")
        success = agent.run(problem, retest_fn=_retest_single_case)
        if success:
            fixed_count += 1

    print(f"\n{'=' * 60}")
    print(f"修复阶段完成: {fixed_count}/{len(selected)} 个问题已修复")
    print(f"{'=' * 60}")


def main():
    print("=" * 60)
    print("MOF 结构组装测试工具")
    print("=" * 60)

    # 1. 扫描参考 CIF
    all_files = scan_reference_cifs()
    if not all_files:
        print("未找到任何参考 CIF 文件，退出。")
        return

    # 2. 交互选择
    selected = select_files(all_files)
    if not selected:
        print("未选中任何文件，退出。")
        return

    print(f"\n将测试 {len(selected)} 个用例。")

    # 3. 逐个处理
    case_details = []  # list of dict: case_name, ref_file, error, gen_path, report

    for i, ref_file in enumerate(selected, 1):
        case_name = os.path.splitext(ref_file)[0]
        print(f"\n{'─' * 60}")
        print(f"[{i}/{len(selected)}] {ref_file}")

        ref_path = os.path.join(os.path.dirname(__file__), "tobacco_output_cifs", ref_file)

        gen_path, error = run_assembly_for_case(ref_file)

        if error:
            print(f"  ⚠️ 跳过: {error}")
            skip_report = StructureReport()
            skip_report.add("组装状态", StructureReport.STATUS_WARN, error)
            case_details.append({
                "case_name": case_name,
                "ref_file": ref_file,
                "error": error,
                "gen_path": None,
                "report": skip_report,
            })
            continue

        # 4. 对比
        report = compare_structures(ref_path, gen_path)
        case_details.append({
            "case_name": case_name,
            "ref_file": ref_file,
            "error": None,
            "gen_path": gen_path,
            "report": report,
        })
        print(report.summary())

    # 5. 汇总报告
    if case_details:
        reports_for_merge = [(cd["case_name"], cd["report"]) for cd in case_details]
        merged, bug_groups = merge_reports(reports_for_merge)
        print("\n" + merged)

        # 保存报告（带时间戳）
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_path = os.path.join(os.path.dirname(__file__), "tmp_output", f"test_report_{ts}.txt")
        os.makedirs(os.path.dirname(report_path), exist_ok=True)
        with open(report_path, "w") as f:
            f.write(merged)
        print(f"\n报告已保存至: {report_path}")

        # 6. 问题分析与修复阶段
        _interactive_fix_loop(case_details, bug_groups)


if __name__ == "__main__":
    main()