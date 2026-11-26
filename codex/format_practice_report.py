from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Mm, Pt


INPUT = Path("/Users/jing/Desktop/专业实践（中期）/专业实践总结报告.docx")
BACKUP = Path("/private/tmp/专业实践总结报告.before_format.docx")

PURPOSE_PARAGRAPHS = [
    "本次专业实践以金属有机框架（Metal–Organic Framework，MOF）对称性感知结构枚举与可视化系统的设计和实现为主线，面向材料计算研究中候选结构生成、筛选和结果管理的实际需求，完成从需求分析、原型搭建到功能实现和应用资料整理的实践。通过将计算方法和软件工程方法结合起来，训练了针对真实科研问题进行需求拆解、模块设计、程序实现、运行调试和成果整理的综合能力。",
    "在系统实现方面，基于 Python 编程语言和 Streamlit（Python Web 应用开发框架）搭建浏览器端操作界面，支持选择 MOF 拓扑、自动分析节点度数需求、配置配体和边缘分子，并调用结构计算流程完成对称性分析、结构枚举、结果验证和晶体学信息文件（Crystallographic Information File，CIF）生成。系统进一步通过 py3Dmol（三维分子结构可视化库）提供候选结构的交互式三维查看，以统计指标、结构汇总表和排名信息辅助结果比较，并支持 CIF 文件和 JavaScript 对象表示法（JavaScript Object Notation，JSON）数据导出。",
    "本次实践的意义在于把相对独立的算法流程组织为可操作、可展示、可归档的 Web 系统，为后续结构分析、第一性原理计算和材料性质研究提供候选结构准备与数据输出支持。同时，实践延续了在航天软件评测中心开展的国产化环境部署、安全配置、MooseFS 分布式文件系统运维与数据分析工作，提升了在国产软硬件环境和关键行业场景下进行系统部署、运行维护和安全管理的综合工程能力，促进了计算机技术与材料科学研究的交叉融合。",
]

PLAN_ROWS = [
    ["2025/3-2025/5", "航天软件评测中心", "航天软件评测中心", "国产化环境部署与安全配置", "已完成"],
    ["2025/6-2025/8", "航天软件评测中心", "航天软件评测中心", "MooseFS 集群运维与数据分析", "已完成"],
    ["2025/9-2025/12", "北京化工大学", "北京化工大学", "完成需求分析与原型搭建", "已完成"],
    ["2026/1-2026/3", "北京化工大学", "北京化工大学", "实现结构枚举与结果展示", "已完成"],
    ["2026/4-2026/6", "北京化工大学", "北京化工大学", "完善三维可视化与数据导出", "已完成"],
]

RESULTS = [
    "1. 完成金属有机框架（MOF）对称性感知结构枚举与可视化系统的 Web 界面开发。基于 Streamlit 组织拓扑选择、配体配置、边缘分子选择、枚举参数设置和结果展示等功能，使用户能够在浏览器中完成候选结构生成、查看和导出的基本流程。",
    "2. 实现拓扑需求自动分析和参数校验功能。系统读取目标拓扑图数据，自动分析不同节点度数的需求，并按度数分组展示可用配体；同时提供边缘分子选择、前 K 个（Top-K）CIF 数量和最大枚举结果数设置，并对缺少配体、拓扑文件不存在等情况给出明确提示。",
    "3. 集成对称性分析、结构枚举、结果验证和结构生成流程。系统调用项目中的枚举与组装模块，输出稳定子大小、轨道大小、有效结构数量等结果信息，并将排名靠前的候选结构生成标准 CIF 文件，为后续结构检查和计算分析提供输入。",
    "4. 完成交互式三维结构展示和结果浏览功能。系统通过 py3Dmol 在浏览器中显示生成的 MOF 晶体结构，支持球棍、棍、球和线框等渲染方式切换；同时以统计指标卡片、有效结构汇总表和排名信息辅助比较不同候选结构的空间形态与结构属性。",
    "5. 完成结构化数据导出及配套资料整理。系统支持逐个下载候选结构的 CIF 文件，并支持预览和下载包含元信息及全部结构详情的 JSON 结果文件；在此基础上整理系统操作手册、源代码材料和软件著作权申请资料，形成较完整的软件交付与使用说明。",
    "6. 作为前期工程实践的补充成果，完成银河麒麟 V10 及配套基础软件的虚拟机模板制作与软件仓库管理，参与部署具备高可用能力的 MooseFS 集群，开发存储输入/输出（Input/Output，I/O）指标采集和日志分析脚本，并在嵌入式测试软件中增加度量功能及配套用户界面（User Interface，UI），进一步积累系统部署、运维分析和测试工具开发经验。",
]

CONCLUSION = (
    "通过本次专业实践，我进一步认识到科研软件开发不仅需要完成算法功能，还需要兼顾需求表达、交互设计、异常处理、结果验证、数据归档和使用说明等工程环节。MOF 结构枚举与可视化系统的实现，使我能够将计算机技术专业知识应用于材料科学研究场景；前期国产化环境和分布式存储实践，则增强了我对系统部署、运维和安全要求的理解。后续将继续围绕系统稳定性、结构验证和面向科研用户的功能完善开展工作。"
)


def delete_paragraph(paragraph) -> None:
    element = paragraph._element
    element.getparent().remove(element)


def set_fonts(run, size: float, *, bold: bool = False, english_only: bool = False) -> None:
    run.font.size = Pt(size)
    run.bold = bold
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts")
        rpr.insert(0, rfonts)
    if english_only:
        rfonts.set(qn("w:ascii"), "Times New Roman")
        rfonts.set(qn("w:hAnsi"), "Times New Roman")
        rfonts.set(qn("w:cs"), "Times New Roman")
    else:
        rfonts.set(qn("w:ascii"), "Times New Roman")
        rfonts.set(qn("w:hAnsi"), "Times New Roman")
        rfonts.set(qn("w:eastAsia"), "宋体")
        rfonts.set(qn("w:cs"), "Times New Roman")


def format_body_paragraph(paragraph, *, first_indent: bool = True) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    pf = paragraph.paragraph_format
    pf.line_spacing = Pt(20)
    pf.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)
    pf.first_line_indent = Mm(7.4) if first_indent else Mm(0)
    for run in paragraph.runs:
        set_fonts(run, 12)


def replace_text(paragraph, text: str, *, size: float = 12, bold: bool = False) -> None:
    if paragraph.runs:
        run = paragraph.runs[0]
        run.text = text
        for extra in paragraph.runs[1:]:
            extra._element.getparent().remove(extra._element)
    else:
        run = paragraph.add_run(text)
    set_fonts(run, size, bold=bold)


def set_cell_margins(cell, top=70, start=80, bottom=70, end=80) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for name, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{name}"))
        if node is None:
            node = OxmlElement(f"w:{name}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def write_cell(cell, text: str, *, center: bool, header: bool = False) -> None:
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER if center else WD_ALIGN_PARAGRAPH.LEFT
    paragraph.paragraph_format.line_spacing = Pt(14)
    paragraph.paragraph_format.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    paragraph.paragraph_format.space_before = Pt(0)
    paragraph.paragraph_format.space_after = Pt(0)
    run = paragraph.add_run(text)
    set_fonts(run, 10.5, bold=header)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
    set_cell_margins(cell)
    if header:
        set_cell_shading(cell, "D9EAF7")


def set_table_borders(table, color="808080", size="6") -> None:
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.first_child_found_in("w:tblBorders")
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        element = borders.find(qn(f"w:{edge}"))
        if element is None:
            element = OxmlElement(f"w:{edge}")
            borders.append(element)
        element.set(qn("w:val"), "single")
        element.set(qn("w:sz"), size)
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), color)


def write_table_title(document, table) -> None:
    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.paragraph_format.space_before = Pt(4)
    title.paragraph_format.space_after = Pt(4)
    title.paragraph_format.line_spacing = Pt(14)
    title.paragraph_format.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    zh = title.add_run("表2-1 实践计划执行情况")
    set_fonts(zh, 10.5, bold=True)
    en = title.add_run("  Table 2-1 Execution of the Professional Practice Plan")
    set_fonts(en, 10.5, bold=True, english_only=True)
    table._tbl.addprevious(title._p)


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)
    shutil.copy2(INPUT, BACKUP)
    document = Document(str(INPUT))
    section = document.sections[0]
    section.page_width = Mm(210)
    section.page_height = Mm(297)
    section.top_margin = Mm(25)
    section.bottom_margin = Mm(25)
    section.left_margin = Mm(27)
    section.right_margin = Mm(27)

    headings = {p.text.strip(): p for p in document.paragraphs if p.text.strip()}
    purpose = next(p for p in document.paragraphs if p.text.startswith("本次专业实践以"))
    heading2 = headings["二、实践计划的执行情况"]
    heading3 = headings["三、实践的主要成果"]
    heading4 = headings["四、专业实践总结报告"]
    table = document.tables[0]

    # Replace the purpose section with three readable, standards-compliant paragraphs.
    replace_text(purpose, PURPOSE_PARAGRAPHS[0])
    format_body_paragraph(purpose)
    for text in PURPOSE_PARAGRAPHS[1:]:
        new_p = heading2.insert_paragraph_before()
        replace_text(new_p, text)
        format_body_paragraph(new_p)

    # Apply heading formatting and make the old placeholder section a real conclusion.
    for heading in (headings["一、本次实践的目的及意义"], heading2, heading3):
        heading.alignment = WD_ALIGN_PARAGRAPH.LEFT
        heading.paragraph_format.space_before = Pt(8)
        heading.paragraph_format.space_after = Pt(4)
        for run in heading.runs:
            set_fonts(run, 14, bold=True)
    replace_text(heading4, "四、实践总结", size=14, bold=True)
    heading4.alignment = WD_ALIGN_PARAGRAPH.LEFT
    heading4.paragraph_format.space_before = Pt(8)
    heading4.paragraph_format.space_after = Pt(4)
    conclusion = heading4.insert_paragraph_before()
    replace_text(conclusion, CONCLUSION)
    format_body_paragraph(conclusion)

    # Remove the old instruction line that followed the placeholder section.
    for p in list(document.paragraphs):
        if p.text.strip() == "要求请参照本专业领域培养方案，可另附页。":
            delete_paragraph(p)

    # Apply the body format to the six results paragraphs and keep their numbering as text.
    for p in document.paragraphs:
        if p.text[:2] in {"1.", "2.", "3.", "4.", "5.", "6."}:
            format_body_paragraph(p, first_indent=False)

    # Make the execution table bilingual and standards-compliant.
    headers = [
        "实践起止时间\nPractice Period",
        "实践基地名称\nPractice Base",
        "地点\nLocation",
        "实践内容\nPractice Content",
        "完成情况\nStatus",
    ]
    for i, cell in enumerate(table.rows[0].cells):
        write_cell(cell, headers[i], center=True, header=True)
    for row in table.rows[1:]:
        for i, cell in enumerate(row.cells):
            write_cell(cell, cell.text, center=i != 3)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    set_table_borders(table)
    write_table_title(document, table)

    # Keep all non-cover content in the specified body format, including the
    # newly inserted purpose and conclusion paragraphs.
    for p in document.paragraphs:
        if any(p._p is q._p for q in (purpose, heading2, heading3, heading4, conclusion)):
            continue
        if p.text.strip() and not p.text.startswith("北京化工大学") and not p.text.startswith("专业学位研究生"):
            if not p.text.startswith("*") and not p.text.startswith("姓") and not p.text.startswith("学") and not p.text.startswith("培养") and not p.text.startswith("专") and not p.text.startswith("校") and not p.text.startswith("填"):
                if p.text.startswith("表2-1"):
                    continue
                if not p.text[:2] in {"1.", "2.", "3.", "4.", "5.", "6."}:
                    format_body_paragraph(p)

    document.save(str(INPUT))
    print(f"Updated and formatted: {INPUT}")
    print(f"Backup: {BACKUP}")


if __name__ == "__main__":
    main()
