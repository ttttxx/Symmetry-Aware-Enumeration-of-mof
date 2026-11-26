from __future__ import annotations

import shutil
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_AUTO_SHAPE_TYPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Pt


SRC = Path("/Users/jing/Desktop/专业实践（中期）/专业实践终期答辩-田婧.pptx")
OUT = Path("/Users/jing/PyCharmMiscProject/mof/codex/专业实践终期答辩-田婧_MOF重排版.pptx")

BLUE = RGBColor(60, 79, 129)
TEXT = RGBColor(48, 48, 48)
MUTED = RGBColor(95, 95, 95)
WHITE = RGBColor(255, 255, 255)


def fit_text(shape, text: str, *, size: float, bold: bool = False, color: RGBColor = TEXT, align=PP_ALIGN.LEFT) -> None:
    text_frame = shape.text_frame
    text_frame.clear()
    text_frame.word_wrap = True
    text_frame.margin_left = Pt(6)
    text_frame.margin_right = Pt(6)
    text_frame.margin_top = Pt(4)
    text_frame.margin_bottom = Pt(4)
    p = text_frame.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = color


def set_body(shape, title: str, body: str, *, title_size: float = 18, body_size: float = 12.5) -> None:
    tf = shape.text_frame
    tf.clear()
    tf.word_wrap = True
    tf.margin_left = Pt(6)
    tf.margin_right = Pt(6)
    tf.margin_top = Pt(4)
    tf.margin_bottom = Pt(4)

    p1 = tf.paragraphs[0]
    r1 = p1.add_run()
    r1.text = title
    r1.font.size = Pt(title_size)
    r1.font.bold = True
    r1.font.color.rgb = TEXT

    p2 = tf.add_paragraph()
    p2.space_before = Pt(5)
    r2 = p2.add_run()
    r2.text = body
    r2.font.size = Pt(body_size)
    r2.font.color.rgb = MUTED


def add_callout(slide, left, top, width, height, title: str, body: str) -> None:
    box = slide.shapes.add_shape(MSO_AUTO_SHAPE_TYPE.ROUNDED_RECTANGLE, left, top, width, height)
    box.fill.solid()
    box.fill.fore_color.rgb = WHITE
    box.line.color.rgb = BLUE
    box.line.width = Pt(1.2)
    tf = box.text_frame
    tf.clear()
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_left = Pt(10)
    tf.margin_right = Pt(10)
    tf.margin_top = Pt(6)
    tf.margin_bottom = Pt(6)

    p1 = tf.paragraphs[0]
    r1 = p1.add_run()
    r1.text = title
    r1.font.size = Pt(14)
    r1.font.bold = True
    r1.font.color.rgb = BLUE

    p2 = tf.add_paragraph()
    p2.space_before = Pt(3)
    r2 = p2.add_run()
    r2.text = body
    r2.font.size = Pt(10.5)
    r2.font.color.rgb = TEXT


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SRC, OUT)
    prs = Presentation(str(OUT))

    # Slide 4: problem statement
    s4 = prs.slides[3]
    fit_text(s4.shapes[1], "问题提出：MOF 候选生成真正卡在两处", size=24, bold=True, color=BLUE)
    fit_text(s4.shapes[8], "组合空间不可直接穷举", size=20, bold=True, color=TEXT, align=PP_ALIGN.CENTER)
    fit_text(
        s4.shapes[4],
        "同一拓扑下，节点配体组合会按指数增长；若直接穷举，低对称构型会迅速淹没搜索，结果里还会夹带大量对称等价重复。",
        size=13,
        color=MUTED,
    )
    fit_text(s4.shapes[15], "图上可行不等于晶体可用", size=20, bold=True, color=TEXT, align=PP_ALIGN.CENTER)
    fit_text(
        s4.shapes[11],
        "MOF 构建单元具有连接几何和取向约束。仅在图上得到着色并不够，还必须把节点和连接体按方向放回周期晶胞，才能得到可审查、可导出的 CIF。",
        size=13,
        color=MUTED,
    )
    fit_text(s4.shapes[6], "Part 01", size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s4.shapes[13], "Part 02", size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Slide 5: why it is hard
    s5 = prs.slides[4]
    fit_text(s5.shapes[0], "难点拆解：从拓扑输入到可复核输出要跨四道坎", size=23, bold=True, color=BLUE)
    fit_text(s5.shapes[4], "拓扑表达", size=18, bold=True, color=TEXT, align=PP_ALIGN.CENTER)
    fit_text(
        s5.shapes[5],
        "周期网络要先约化为 primitive cell 上的 quotient graph，并把 topology、节点度数和配体池组织成可计算的着色空间。",
        size=12.2,
        color=MUTED,
    )
    fit_text(s5.shapes[10], "无冗余枚举", size=18, bold=True, color=TEXT, align=PP_ALIGN.CENTER)
    fit_text(
        s5.shapes[11],
        "要在真旋转群作用下只输出 symmetry-distinct 着色；子群格下降、Canonical DFS 与按对称性优先的调度缺一不可。",
        size=12.2,
        color=MUTED,
    )
    fit_text(s5.shapes[16], "晶体重构", size=18, bold=True, color=TEXT, align=PP_ALIGN.CENTER)
    fit_text(
        s5.shapes[17],
        "抽象着色还要经过 Kabsch 刚体对齐、连接体轴向旋转与坐标归一化，否则结果无法真正落到晶体模型。",
        size=12.2,
        color=MUTED,
    )
    fit_text(s5.shapes[19], "可信交付", size=18, bold=True, color=TEXT, align=PP_ALIGN.CENTER)
    fit_text(s5.shapes[22], "3D 视图、结构表、稳定子/轨道属性、CIF/JSON 与运行日志必须放在同一链路里，结果才可复核。", size=12.5, color=MUTED)
    fit_text(s5.shapes[24], "不能只给“结构数”，还要给验证信息和后续计算可直接使用的文件。", size=12.5, color=MUTED)
    fit_text(s5.shapes[26], "研究人员最终需要的是“能继续算、能回头查、能解释为什么保留这批候选”的结果。", size=12.5, color=TEXT)

    # Slide 6: technical route
    s6 = prs.slides[5]
    fit_text(s6.shapes[0], "解决思路：四个模块协同闭环", size=24, bold=True, color=BLUE)
    fit_text(
        s6.shapes[21],
        "我把论文里的“枚举—重构—验证”流程真正做成了可操作的软件系统，而不是只停留在离线脚本。",
        size=13,
        color=TEXT,
    )
    fit_text(s6.shapes[5], "读取拓扑图后计算自同构群，只保留与后续重构兼容的真旋转置换。", size=12.2, color=MUTED)
    fit_text(s6.shapes[6], "对称性分析", size=17, bold=True, color=TEXT)
    fit_text(s6.shapes[7], "01", size=20, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s6.shapes[8], "按轨道数从高对称到低对称调度任务，Canonical DFS 在搜索中直接剪枝去重。", size=12.2, color=MUTED)
    fit_text(s6.shapes[9], "子群格枚举", size=17, bold=True, color=TEXT)
    fit_text(s6.shapes[10], "02", size=20, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s6.shapes[11], "用 stabilizer / orbit 统计和 Burnside 恒等式检查覆盖率，并与 CASM 结果交叉对账。", size=12.2, color=MUTED)
    fit_text(s6.shapes[12], "正确性核对", size=17, bold=True, color=TEXT)
    fit_text(s6.shapes[13], "03", size=20, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s6.shapes[14], "把结构表、运行日志、3D 视图与 CIF/JSON 下载接到同一前端，支持直接复核。", size=12.2, color=MUTED)
    fit_text(s6.shapes[15], "界面交付", size=17, bold=True, color=TEXT)
    fit_text(s6.shapes[16], "04", size=20, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Slide 7: deep technical explanation over the workflow figure
    s7 = prs.slides[6]
    fit_text(s7.shapes[0], "核心方法：从商图着色到晶体重构", size=24, bold=True, color=BLUE)
    fit_text(s7.shapes[2], "商图约化 + 子群格下降 + Kabsch 刚体重构", size=18, bold=True, color=BLUE, align=PP_ALIGN.CENTER)
    add_callout(
        s7,
        720000,
        1180000,
        2650000,
        1130000,
        "1 商图约化",
        "先把周期平移约掉，在 primitive cell 上把 MOF 设计问题改写成 building block 着色问题。",
    )
    add_callout(
        s7,
        8870000,
        1180000,
        2950000,
        1230000,
        "2 对称性指标",
        "用 rho = stabilizer_size / |G| 衡量对称层级；轨道越小，说明构型越高对称，优先输出。",
    )
    add_callout(
        s7,
        780000,
        4880000,
        3100000,
        1320000,
        "3 枚举策略",
        "任务按 全群 → 真子群 → 恒等层 调度；必要时只截断 identity tail，保留完整高对称前缀。",
    )
    add_callout(
        s7,
        8620000,
        4800000,
        3250000,
        1420000,
        "4 晶体重构",
        "节点用 Kabsch 对齐连接位点，连接体沿边向量放置并做 RMSD 微调，最终生成 sanitizable CIF。",
    )

    # Slide 8: screenshot chain as result
    s8 = prs.slides[7]
    fit_text(s8.shapes[0], "成果一：系统已具备完整功能链", size=23, bold=True, color=BLUE)
    fit_text(s8.shapes[17], "真实界面截图，展示一次运行从参数输入到结果归档的完整流程。", size=13, color=TEXT, align=PP_ALIGN.CENTER)
    fit_text(s8.shapes[3], "参数配置", size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s8.shapes[10], "选择拓扑、节点配体、边缘分子以及 Top-K / max results 等关键参数。", size=12, color=MUTED)
    fit_text(s8.shapes[11], "动态需求", size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s8.shapes[12], "按节点度数自动生成可选配体池，避免化学计量约束和前端输入脱节。", size=12, color=MUTED)
    fit_text(s8.shapes[13], "三维审查", size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s8.shapes[14], "Top-K 结构可并行渲染，并支持 Ball-and-Stick / Stick / Sphere / Line 样式切换。", size=12, color=MUTED)
    fit_text(s8.shapes[15], "结果归档", size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s8.shapes[16], "候选结构可逐个导出 CIF，完整运行结果可导出 JSON，便于复算与论文留档。", size=12, color=MUTED)

    # Slide 9: main interface
    s9 = prs.slides[8]
    fit_text(s9.shapes[0], "成果二：主界面把算法指标和结构结果放在同一链路上", size=22, bold=True, color=BLUE)
    fit_text(s9.shapes[2], "系统主界面", size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s9.shapes[3], "参数—计算—结果", size=14, bold=True, color=TEXT, align=PP_ALIGN.CENTER)
    fit_text(s9.shapes[6], "动态拓扑需求", size=15, bold=True, color=TEXT)
    fit_text(s9.shapes[8], "读取 topology requirement，并按 3-c / 4-c 等槽位组织节点配体选择。", size=11.5, color=MUTED)
    fit_text(s9.shapes[12], "枚举与验证", size=15, bold=True, color=TEXT)
    fit_text(s9.shapes[14], "前端可直接看到 |G|、Valid Structures、CIF 数量和运行日志，而不是只拿到黑箱输出。", size=11.5, color=MUTED)
    fit_text(s9.shapes[18], "三维结果", size=15, bold=True, color=TEXT)
    fit_text(s9.shapes[20], "候选结构可旋转、缩放并切换渲染样式，便于快速排查方向性和局部连接异常。", size=11.5, color=MUTED)
    fit_text(s9.shapes[24], "结构属性", size=15, bold=True, color=TEXT)
    fit_text(s9.shapes[26], "表格同步给出 coloring、stabilizer size 和 orbit size，结果能按对称性解释。", size=11.5, color=MUTED)
    fit_text(s9.shapes[30], "完整归档", size=15, bold=True, color=TEXT)
    fit_text(s9.shapes[32], "JSON 输出保留参数、验证信息和全部结构属性，后续计算和论文整理都能复现。", size=11.5, color=MUTED)
    fit_text(s9.shapes[36], "工程边界", size=15, bold=True, color=TEXT)
    fit_text(s9.shapes[38], "系统面向候选结构准备和复核，不直接替代 DFT 或独立稳定性计算。", size=11.5, color=MUTED)

    # Slide 10: final results summary with technical weight
    s10 = prs.slides[9]
    fit_text(s10.shapes[0], "成果三：终期交付既有算法内核，也有可复核证据", size=23, bold=True, color=BLUE)
    fit_text(s10.shapes[1], "可穷举规模", size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s10.shapes[2], "约束空间可处理到约 5×10^5 个配置，并保持按对称性优先输出。", size=11.5, color=MUTED)
    fit_text(s10.shapes[7], "正确性核对", size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s10.shapes[8], "穷举案例与 Burnside / CASM 对账一致；受限运行时也能报告覆盖率边界。", size=11.5, color=MUTED)
    fit_text(s10.shapes[13], "三维审查", size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s10.shapes[14], "Top-K 候选结构支持并行查看，便于从图级结果过渡到晶体级检查。", size=11.5, color=MUTED)
    fit_text(s10.shapes[16], "结构导出", size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s10.shapes[17], "接受的着色可恢复为 periodic CIF，单结构文件可直接下载用于后续计算。", size=11.5, color=MUTED)
    fit_text(s10.shapes[10], "运行留痕", size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s10.shapes[11], "计算过程和验证日志可回看，避免一次运行结束后只剩结果文件、无法追溯。", size=11.5, color=MUTED)
    fit_text(s10.shapes[4], "数据归档", size=15, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    fit_text(s10.shapes[5], "完整 JSON 保留元数据、验证信息和结构属性，保证一次运行可复算、可留档。", size=11.5, color=MUTED)
    add_callout(
        s10,
        3980000,
        840000,
        4300000,
        700000,
        "总结",
        "终期成果不是“做了一个界面”，而是把对称枚举、晶体重构、正确性核对和科研交付串成了闭环。",
    )

    prs.save(str(OUT))
    print(OUT)


if __name__ == "__main__":
    main()
