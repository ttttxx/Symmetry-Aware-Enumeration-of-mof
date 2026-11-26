from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Mm, Pt


SOURCE = Path("/Users/jing/Desktop/专业实践（中期）/专业实践总结报告.docx")
OUTPUT = Path("/Users/jing/PyCharmMiscProject/mof/codex/专业实践总结报告_第四部分重写_2026-07-13.docx")
ASSET_DIR = Path("/Users/jing/PyCharmMiscProject/mof/codex/practice_report_assets")


BLOCKS = [
    ("h1", "1、专业实践基本情况"),
    ("h2", "1.1 实践背景"),
    (
        "p",
        "金属有机框架材料是一类由金属节点与有机配体组装形成的晶态多孔材料。它的吸引力并不只在于结构好看，"
        "更在于“可设计”三个字：节点可以换，配体可以换，拓扑也可以换，孔道大小、化学环境和潜在功能随之发"
        "生明显变化。正因为这种可调性，MOF 在气体储存与分离、催化、传感和分子识别等方向持续受到关注。与此"
        "同时，材料设计面对的不是几个离散候选，而是一个规模极大的组合空间。只要候选构件库稍微扩大，可能的"
        "结构数目就会迅速增长。研究工作因此很容易卡在第一步：不是没有思路，而是候选太多，人工试配既慢也不"
        "可靠。"
    ),
    (
        "p",
        "从计算角度看，MOF 候选结构生成至少有两道难关。第一道难关是组合爆炸。相同拓扑下，不同构件在不同位"
        "点上的放置会产生海量方案，其中大量方案只是彼此在对称操作下等价，若不加区分地全部生成，计算会被冗"
        "余结构淹没。第二道难关是“图上成立”并不自动等于“结构上可用”。有些方案在抽象连接关系层面没有问题，"
        "但一旦落实到真实构件的朝向、连接点匹配和晶胞尺度，几何冲突就会出现。也就是说，候选结构生成不能只"
        "停留在图论层面的排列组合，还要考虑后续晶体重建和结构检查能否顺利接上。"
    ),
    (
        "p",
        "本次在北京化工大学开展的专业实践，就是围绕这一研究痛点展开的。课题组前期已经形成了拓扑读取、对称"
        "性分析、结构枚举、验证和晶体重建等脚本，但这些功能主要面向命令行实验，适合熟悉代码的使用者，不适"
        "合高频迭代的科研场景。研究人员往往需要反复切换脚本、结果文件和外部可视化工具，参数稍一变化，就要"
        "重新组织一遍流程。随着课题推进，系统化整合变得很有必要：需要一个能够直接面向研究过程的 Web 界面，"
        "把拓扑选择、构件配置、枚举计算、结果查看、结构预览和文件导出放进同一套工作流里。"
    ),
    (
        "p",
        "需要说明的是，这次专业实践并不是从零开始的一段孤立开发。在前期航天软件评测中心的实践中，我已经参"
        "与过国产化环境部署、安全配置、基础软件安装、虚拟机模板制作、MooseFS 分布式存储运维和日志分析脚本"
        "开发等工作。那一阶段的内容与 MOF 研究不直接相连，但它让我对软件环境一致性、运行稳定性、错误定位"
        "和工程交付有了更具体的认识。正是有了这部分积累，后续在北京化工大学做科研软件整合时，我才会更自然"
        "地把异常提示、日志留存、输入约束和结果可复用性放到与功能实现同等重要的位置。"
    ),
    ("h2", "1.2 实践任务与目标"),
    (
        "p",
        "围绕北京化工大学阶段的主要工作，本次实践的核心任务可以概括为一句话：在不改变项目科学边界的前提下，"
        "把已有的 MOF 对称性感知结构枚举与重建流程做成一个可直接使用、可直接展示、也便于后续整理资料的 Web"
        " 系统。为完成这一任务，我需要同时处理三类问题。其一是交互问题，即如何把研究人员熟悉的操作步骤转化"
        "为清晰的页面逻辑；其二是算法接入问题，即界面调用的必须是真实的对称性分析、枚举、验证和重建模块，"
        "不能为了好展示而把流程简化成示意版；其三是交付问题，即运行后的结果不仅要看得见，还要能下载、能复"
        "用、能写进操作手册和软件著作权申请材料。"
    ),
    (
        "p",
        "据此，我将技术路线确定为“拓扑数据读取与需求分析—对称性提取—不等价构型枚举—几何重建与验证—结果"
        "展示与导出”五个连续环节。界面层由 Streamlit 组织，流程控制由 main.py 负责；symmetry.py 从拓扑图中提"
        "取真旋转诱导的置换群，enumerator.py 在对子群格分层的基础上完成对称性不等价着色搜索，assembler.py"
        " 负责把抽象构型重建为可保存的 CIF 晶体结构，最终再由 py3Dmol 和结果表格把生成结果反馈给用户。这条"
        "路线的关键，不是把很多模块硬拼到一起，而是保证每个环节的输入输出能够自然衔接。"
    ),
    ("h1", "2、专业实践过程与实施内容"),
    ("h2", "2.1 前期工程实践过程"),
    (
        "p",
        "前期在航天软件评测中心的实践，主要工作并不在材料软件本身，而在系统基础环境。实践中，我参与了银河"
        "麒麟 V10 环境部署、安全配置、JDK 与 Python 等基础软件安装、虚拟机模板制作以及软件仓库管理工作。工"
        "作内容看起来偏底层，但它直接影响后续软件能否稳定运行。这个阶段让我更清楚地认识到，工程系统的问题"
        "往往不是“有没有代码”，而是“代码放到具体环境里以后能不能连续跑、别人能不能照着复现、出了问题能不"
        "能快速定位”。"
    ),
    (
        "p",
        "除此之外，我还参与了 MooseFS 分布式文件系统集群的部署与高可用维护，编写了存储 I/O 指标采集和日志分"
        "析脚本，并在嵌入式测试软件中补充了度量功能和界面。这些工作训练了我对资源瓶颈、运行状态和异常信息"
        "的敏感性。后来再看 MOF Web 系统的开发，会发现很多处理思路并不是凭空出现的，例如为什么要做前置输入"
        "检查、为什么要保留运行日志、为什么要对结果数做上限控制，本质上都与前期形成的工程意识有关。"
    ),
    ("h2", "2.2 MOF系统需求分析与交互设计"),
    (
        "p",
        "进入北京化工大学阶段后，我没有直接从页面美化入手，而是先梳理研究人员一次完整使用流程。实际操作通"
        "常是这样的：先选拓扑，再确认当前拓扑需要哪些配位数的节点构件，然后选择候选节点与边缘构件，设定枚"
        "举结果规模，运行计算，最后查看统计指标、结构预览和导出文件。问题恰恰出在这里。原有脚本虽然能跑，"
        "但整个过程分散，用户必须自己记住拓扑约束、自己组织输入，还要在不同文件之间查找结果。于是，我把页"
        "面交互设计成侧边栏负责参数输入，主区域负责状态展示与结果浏览，让“配置—运行—查看”这条链尽量短。"
    ),
    (
        "p",
        "在具体实现上，app.py 会先读取当前拓扑对应的图数据，再自动分析节点度数需求，按配位数动态给出可选构"
        "件。这样做的好处很直接：用户不再需要提前记忆某个拓扑究竟需要几个三配位节点、几个四配位节点，界面"
        "本身就能把约束信息解释清楚。与此同时，我把常见错误尽量前置处理。没有为某一类配位位点选择构件，界"
        "面会立即提示；拓扑文件缺失，程序会明确报错；运行按钮触发后，标准输出会被重定向到日志缓冲区，页面"
        "同步显示“正在运行”的状态和最终信息。这样的页面并不追求花哨，而是尽量贴近真实科研场景中“少切换、"
        "少猜测、少返工”的需求。"
    ),
    ("h2", "2.3 对称性感知结构枚举算法集成"),
    (
        "p",
        "MOF 结构枚举的第一步，是把无限周期晶体网络转化为可计算的有限对象。论文中采用的做法，是把晶体网络"
        "视为无限周期图，再对平移群进行商运算，得到单胞尺度上的有限商图："
    ),
    ("eq", "G = G_periodic / T = (V_cell, E_cell)"),
    ("eq_no", "（4-1）"),
    (
        "p",
        "在代码里，这一步主要对应 main.py 对图数据的读取和增强。程序不仅装载节点、边和晶格信息，还进一步计"
        "算后续重建所需的 alpha_matrix、Bstar 等拓扑辅助量，并统计不同配位数节点的需求数量。于是，原本依赖"
        "人工理解的拓扑约束被转化成了结构化输入，页面可以据此自动生成参数项，枚举器也可以据此进行约束搜索。"
    ),
    (
        "p",
        "在得到有限商图之后，候选结构生成可以表述为节点着色问题。若候选构件集合记为 C，则全部候选配置空间可"
        "写为："
    ),
    ("eq", "X = C^(V_cell)"),
    ("eq_no", "（4-2）"),
    (
        "p",
        "问题在于，式（4-2）给出的空间通常远大于实际需要考察的结果集合。大量候选只是彼此在对称操作下等价。"
        "如果不先消除这部分冗余，后续重建与筛选的成本会迅速上升。因此，系统在真正枚举之前，先由 symmetry.py"
        " 对拓扑结构做对称性分析，提取由真旋转操作诱导的置换群 Γ，而不是直接采用全部组合自同构。这样处理的"
        "目的，是保证后续区分仍然与真实构件的取向重建保持一致。"
    ),
    (
        "p",
        "在群 Γ 的作用下，同一候选构型对应一个轨道，轨道大小与稳定子大小满足轨道—稳定子关系："
    ),
    ("eq", "|O_c| = |Γ| / |Stab(c)|"),
    ("eq_no", "（4-3）"),
    (
        "p",
        "为了便于比较不同构型的对称程度，论文进一步引入了归一化对称性指标："
    ),
    ("eq", "ρ(c) = |Stab(c)| / |Γ| = 1 / |O_c|"),
    ("eq_no", "（4-4）"),
    (
        "p",
        "enumerator.py 的实现并不是对所有节点盲目深搜，而是先按对子群格分层组织搜索任务，优先处理高对称子空"
        "间，再逐步进入低对称区域。每个分层内部，程序为同一轨道上的节点统一赋值；一旦发现某个部分着色在字"
        "典序意义下已经不是当前等价类的最小代表元，就立即剪枝，不再继续扩展。最后，完整着色会被映射到全群"
        "下的规范表示，并插入 Trie 结构中去重。这样得到的不是“全部排列”，而是对称性不等价的候选结构集合。"
    ),
    (
        "p",
        "对于可穷举空间，系统还使用 Burnside 恒等式检查轨道计数的一致性："
    ),
    ("eq", "Σ_{c∈Ω} 1 / |Stab(c)| = |X| / |Γ|"),
    ("eq_no", "（4-5）"),
    (
        "p",
        "式（4-5）左侧由程序对已经得到的不等价结果逐项累加，右侧由配置空间规模和群阶共同确定。当两者一致时，"
        "说明当前约束下的枚举在轨道意义上是完备的。若恒等子群对应的低对称空间过大，代码会启动预算机制，优"
        "先保留高对称部分和更有信息量的前缀结果。这一策略既照顾了理论完整性，也照顾了实际运行成本。"
    ),
    ("h2", "2.4 几何重建、系统集成与结果导出"),
    (
        "p",
        "仅仅得到抽象着色，并不能直接说明一个 MOF 候选结构已经可以使用。真正能够进入后续结构检查或计算的，"
        "必须是带有具体原子坐标、构件朝向和晶格信息的晶体模型。因此，系统在枚举之后还要完成几何重建。"
        "assembler.py 负责这一过程：它根据节点模板、边缘构件、拓扑连接关系和晶格参数，将抽象的“节点—构件”"
        "映射恢复为实际坐标。这里最关键的一个步骤，是采用刚体对齐思想，使模板构件的连接位点尽量贴合目标位"
        "置，其最小二乘目标可写为："
    ),
    ("eq", "min_{R,t} Σ_{i=1}^m ||R x_i + t - y_i||_2^2,  s.t.  R^T R = I,  det(R) = 1"),
    ("eq_no", "（4-6）"),
    (
        "p",
        "式（4-6）对应的正是 Kabsch 对齐问题，其中 x_i 与 y_i 分别表示模板构件和目标连接位点的坐标集合，R 与 t"
        " 表示旋转和平移。围绕这一目标，代码进一步完成了晶胞尺度自适应放缩、节点与边缘连接点融合、桥接键修"
        "正、空间群识别和结构校验。经过这些处理后，系统输出的才是可保存、可预览、也可继续用于后续计算的 CIF"
        " 文件。也正因为有这一阶段，整个系统并不只是一个图论枚举器，而是一套能够落到实际结构文件上的研究工"
        "具。"
    ),
    (
        "p",
        "在系统层面，我把这一后端流程完整接入到了 Web 页面。页面顶部用指标卡给出当前拓扑、自同构群阶数、有"
        "效结构数量和已生成 CIF 数量；中间用表格列出每个候选结构的着色、稳定子大小和轨道大小；对排名靠前的"
        "结果，则用 py3Dmol 提供交互式三维预览，并支持球棍、棍、球和线框等多种渲染方式。用户如果只想快速检"
        "查结构，可以停留在可视化页面；如果要继续分析，也可以直接下载 CIF 或完整 JSON。输出结果因此同时兼顾"
        "了“看得见”和“接得上”。"
    ),
    (
        "p",
        "为了让系统真正可复用，我在工程实现中还加入了若干保障机制。例如，缺少必要输入时在运行前阻断，避免"
        "浪费计算；对最大结果数与 Top-K CIF 数量设置控制参数，避免界面因结果过多而失去可读性；对运行日志做"
        "集中展示，便于定位问题；对异常堆栈保留页面回显，减少调试成本。这些处理不直接增加论文里的公式，却直"
        "接决定了系统在真实研究场景中是否耐用。"
    ),
    ("table", "表4-1 北京化工大学阶段 MOF 系统主要模块与实践内容"),
    ("h1", "3、专业实践成果"),
    ("h2", "3.1 系统界面与交互实现效果"),
    (
        "p",
        "本次实践最直接、也最容易被验证的成果，是完成了 MOF 对称性感知结构枚举与可视化系统的 Web 界面。与"
        "原来分散在脚本、结果文件和外部查看工具中的流程相比，新的系统把参数配置、运行触发和结果浏览组织到"
        "了同一页面中。用户首先在左侧选择拓扑，系统会自动给出当前拓扑的节点配位需求；随后根据配位数动态展"
        "开可选构件，用户无需预先记忆约束条件，就能直接完成输入配置。"
    ),
    (
        "p",
        "图4-1 展示了系统主界面在运行过程中的状态。可以看到，左侧侧边栏保留了拓扑选择、节点构件配置、边缘"
        "构件指定、结果规模控制等核心输入项；右侧主区域则展示运行提示与日志。这样的页面组织不是为了“好看”，"
        "而是为了减少科研人员在命令行、日志文件和结果目录之间来回切换的成本。图4-2 则进一步显示了动态构件"
        "选择的细节：当前拓扑需要的 3-c 与 4-c 节点数目会被自动识别，界面只展示满足该配位数要求的构件候选。"
        "这说明系统已经把底层拓扑约束真正转化为了可操作的交互逻辑。"
    ),
    ("figure", "图4-1 系统主界面与运行日志展示|image11.png|155"),
    ("figure", "图4-2 基于拓扑需求的动态构件选择界面|image10.png|95"),
    ("h2", "3.2 结构生成结果展示与数据输出效果"),
    (
        "p",
        "在结果展示方面，系统并没有停留在“输出一个文件路径”的水平，而是进一步把结构结果做成了可直接浏览的"
        "界面。枚举完成后，页面会给出当前拓扑名称、自同构群阶数、有效结构数与已生成 CIF 数量，并将候选结构"
        "按稳定子大小等指标进行展示。对使用者来说，这些信息的价值很明确：一方面可以迅速了解当前计算规模，"
        "另一方面也能据此判断哪些结构具有更高对称性、哪些结果值得优先查看。"
    ),
    (
        "p",
        "图4-3 展示了系统的三维结构查看页面。前几名候选结构会以并列形式显示，用户可以直接切换渲染方式，对"
        "结构骨架、节点位置和连接形态进行直观检查。对于以往必须借助外部软件才能完成的这一步，如今在浏览器"
        "中就可以完成，明显提高了查看效率。更重要的是，这种三维展示并非脱离算法结果的示意图，而是基于枚举"
        "和重建流程生成的真实 CIF 结构文件。"
    ),
    ("figure", "图4-3 候选 MOF 结构的三维查看与渲染方式切换|image12.png|155"),
    (
        "p",
        "除了三维预览，系统还对全部有效结构给出了表格化结果，并提供 JSON 导出。图4-4 所示界面中，可以直接查"
        "看到每个候选结构的着色方案、稳定子大小和轨道大小。这样的结果组织形式对后续工作非常重要：若只保留"
        "结构文件，很难快速比较不同候选之间的关系；而把统计属性与下载入口统一放在同一页面中，研究人员就可"
        "以在浏览、筛选和归档之间自然切换。JSON 文件则进一步保存了拓扑、验证信息和结构明细，为后续脚本分析、"
        "论文表述和计算任务衔接提供了统一数据源。"
    ),
    ("figure", "图4-4 全部有效结构统计表与 JSON 输出界面|image8.png|160"),
    (
        "p",
        "从实现效果看，系统已经形成了一条较完整的“输入—计算—预览—导出”闭环。它既能服务于单次结果查看，也"
        "能为后续结构检查、几何分析和第一性原理计算准备输入材料。相比仅仅生成若干中间文件，这样的成果更具"
        "备研究支撑价值，也更符合专业实践对于“面向真实场景解决问题”的要求。"
    ),
    ("h2", "3.3 配套资料与交付成果"),
    (
        "p",
        "围绕系统开发，我还同步完成了操作手册、软件著作权申请资料、代码材料和相关说明文本的整理工作。这样"
        "做的意义在于，实践成果不再只是停留在个人电脑上的一组代码，而是形成了“系统实现—使用说明—资料归档”"
        "相互支撑的成果集合。无论是课题组内部复用，还是后续答辩、归档和成果申报，这部分材料都提供了明确依"
        "据。"
    ),
    (
        "p",
        "更具体地说，本次成果至少体现为三个落地点。第一，形成了可运行、可演示、可下载结果的 Web 系统；第二，"
        "形成了带有 CIF 与 JSON 标准化输出的数据组织方式，能够继续服务于后续研究；第三，形成了操作手册和软"
        "著申请资料，使系统具备较完整的交付表达。就专业实践而言，这三点比单纯完成一个页面或单纯写出一段算"
        "法更有分量，因为它们共同说明该工作已经从“代码实现”走向了“可使用的研究软件成果”。"
    ),
    ("h1", "4、实践收获"),
    (
        "p",
        "通过这段实践，我最大的收获不是单一技术点，而是对“科研软件开发”这件事有了更完整的理解。以前更容易"
        "把注意力集中在算法能不能跑通，现在会自然地追问：输入是否清楚、边界是否明确、错误能否定位、结果能"
        "否复用。前期在航天软件评测中心形成的工程意识，与后期 MOF Web 系统开发中的算法整合经验，在这里真正"
        "接到了一起。"
    ),
    (
        "p",
        "另一方面，这次实践也让我更适应跨学科工作方式。材料研究中的拓扑、配位数、空间群、CIF 等概念，必须"
        "先被理解，再被翻译成程序中的数据结构和流程控制，最后才可能变成可用界面。这个过程并不轻松，但它让"
        "我意识到，计算机技术在科研场景中的价值，往往就体现在这种“把复杂问题变成可操作系统”的能力上。"
    ),
    ("h1", "5、意见建议"),
    (
        "p",
        "结合本次实践，我认为自己后续最需要继续加强的，一是对材料领域问题的理解深度，二是对研究软件长期维"
        "护方式的把握。前者决定系统是否真正贴合研究需求，后者决定系统能否在项目推进过程中持续可用。尤其是"
        "面对跨学科任务时，若只关注代码实现，很容易把系统做成“能演示但不好用”的原型。"
    ),
    (
        "p",
        "因此，今后无论继续完善该系统，还是承担类似课题，我都需要把需求梳理、异常处理、结果归档和文档表达"
        "放在更靠前的位置。对自己而言，这不是额外负担，而是把技术工作真正做实所必须补齐的一部分。"
    ),
]


TABLE_ROWS = [
    [
        "拓扑输入与需求分析",
        "app.py、main.py",
        "读取拓扑图数据，自动统计不同配位数节点需求，动态生成构件选择项。",
        "将领域约束直接转化为可操作输入，减少人工配置错误。",
    ],
    [
        "对称性分析与结构枚举",
        "symmetry.py、enumerator.py",
        "提取真旋转诱导的置换群，在对子群格分层基础上执行规范性剪枝搜索。",
        "输出对称性不等价构型，降低冗余计算量。",
    ],
    [
        "几何重建与结构校验",
        "assembler.py",
        "完成刚体对齐、晶胞放缩、结构修正、空间群识别与 CIF 生成。",
        "把抽象枚举结果落实为可继续用于研究的真实结构文件。",
    ],
    [
        "结果展示与数据导出",
        "app.py、py3Dmol",
        "展示统计指标、三维结构、结果表格，并支持 CIF 与 JSON 下载。",
        "提升结果可读性、可复用性和后续流程衔接能力。",
    ],
]


def delete_paragraph(paragraph) -> None:
    element = paragraph._element
    element.getparent().remove(element)


def set_fonts(run, size: float, *, bold: bool = False, east_asia: str = "宋体", ascii_font: str = "Times New Roman") -> None:
    run.font.size = Pt(size)
    run.bold = bold
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts")
        rpr.insert(0, rfonts)
    rfonts.set(qn("w:ascii"), ascii_font)
    rfonts.set(qn("w:hAnsi"), ascii_font)
    rfonts.set(qn("w:eastAsia"), east_asia)
    rfonts.set(qn("w:cs"), ascii_font)


def format_body(paragraph, *, first_indent: bool = True, align=WD_ALIGN_PARAGRAPH.JUSTIFY) -> None:
    paragraph.alignment = align
    pf = paragraph.paragraph_format
    pf.line_spacing = Pt(20)
    pf.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)
    pf.first_line_indent = Mm(7.4) if first_indent else Mm(0)
    for run in paragraph.runs:
        set_fonts(run, 12)


def add_body(document: Document, text: str) -> None:
    p = document.add_paragraph()
    p.add_run(text)
    format_body(p)


def add_heading(document: Document, text: str, level: int) -> None:
    p = document.add_paragraph()
    run = p.add_run(text)
    size = 14 if level == 1 else 12
    set_fonts(run, size, bold=True)
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    pf = p.paragraph_format
    pf.line_spacing = Pt(20)
    pf.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    pf.space_before = Pt(8 if level == 1 else 4)
    pf.space_after = Pt(2)
    pf.first_line_indent = Mm(0)


def set_cell_margins(cell, top=60, start=80, bottom=60, end=80) -> None:
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


def remove_table_borders(table) -> None:
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
        element.set(qn("w:val"), "nil")


def add_equation(document: Document, expr: str, number: str) -> None:
    table = document.add_table(rows=1, cols=2)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    remove_table_borders(table)
    left, right = table.rows[0].cells
    left.width = Mm(126)
    right.width = Mm(18)
    set_cell_margins(left, top=10, bottom=10)
    set_cell_margins(right, top=10, bottom=10)

    p_left = left.paragraphs[0]
    p_left.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pf_left = p_left.paragraph_format
    pf_left.line_spacing = Pt(18)
    pf_left.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    run_left = p_left.add_run(expr)
    set_fonts(run_left, 11.5, east_asia="Times New Roman", ascii_font="Times New Roman")

    p_right = right.paragraphs[0]
    p_right.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    pf_right = p_right.paragraph_format
    pf_right.line_spacing = Pt(18)
    pf_right.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    run_right = p_right.add_run(number)
    set_fonts(run_right, 11, east_asia="宋体")


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


def shade_cell(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def add_module_table(document: Document, title_text: str) -> None:
    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pf = title.paragraph_format
    pf.space_before = Pt(6)
    pf.space_after = Pt(4)
    pf.line_spacing = Pt(16)
    pf.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    run = title.add_run(title_text)
    set_fonts(run, 10.5, bold=True)

    table = document.add_table(rows=1 + len(TABLE_ROWS), cols=4)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    set_table_borders(table)

    headers = ["模块", "对应代码", "主要实现内容", "实践价值"]
    widths_mm = [24, 33, 52, 41]
    for idx, cell in enumerate(table.rows[0].cells):
        cell.width = Mm(widths_mm[idx])
        set_cell_margins(cell, top=80, bottom=80)
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        shade_cell(cell, "D9EAF7")
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        pf = p.paragraph_format
        pf.line_spacing = Pt(14)
        pf.line_spacing_rule = WD_LINE_SPACING.EXACTLY
        run = p.add_run(headers[idx])
        set_fonts(run, 10.5, bold=True)

    for row_idx, values in enumerate(TABLE_ROWS, start=1):
        row = table.rows[row_idx]
        for col_idx, value in enumerate(values):
            cell = row.cells[col_idx]
            cell.width = Mm(widths_mm[col_idx])
            set_cell_margins(cell, top=70, bottom=70)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if col_idx < 2 else WD_ALIGN_PARAGRAPH.LEFT
            pf = p.paragraph_format
            pf.line_spacing = Pt(14)
            pf.line_spacing_rule = WD_LINE_SPACING.EXACTLY
            run = p.add_run(value)
            set_fonts(run, 10)


def add_figure(document: Document, spec: str) -> None:
    caption, filename, width_mm = spec.split("|")
    image_path = ASSET_DIR / filename
    if not image_path.exists():
        raise FileNotFoundError(image_path)

    p_img = document.add_paragraph()
    p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pf_img = p_img.paragraph_format
    pf_img.space_before = Pt(4)
    pf_img.space_after = Pt(2)
    pf_img.line_spacing = Pt(16)
    pf_img.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    run = p_img.add_run()
    run.add_picture(str(image_path), width=Mm(float(width_mm)))

    p_cap = document.add_paragraph()
    p_cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pf_cap = p_cap.paragraph_format
    pf_cap.space_before = Pt(0)
    pf_cap.space_after = Pt(4)
    pf_cap.line_spacing = Pt(16)
    pf_cap.line_spacing_rule = WD_LINE_SPACING.EXACTLY
    run_cap = p_cap.add_run(caption)
    set_fonts(run_cap, 10.5)


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SOURCE, OUTPUT)
    document = Document(str(OUTPUT))

    section = document.sections[0]
    section.page_width = Mm(210)
    section.page_height = Mm(297)
    section.top_margin = Mm(25)
    section.bottom_margin = Mm(25)
    section.left_margin = Mm(27)
    section.right_margin = Mm(27)

    heading4 = None
    for paragraph in document.paragraphs:
        if paragraph.text.strip() == "四、专业实践总结报告":
            heading4 = paragraph
            break
    if heading4 is None:
        raise ValueError("未找到“四、专业实践总结报告”标题。")

    paragraphs = list(document.paragraphs)
    start_idx = next(i for i, p in enumerate(paragraphs) if p._p is heading4._p)
    for paragraph in paragraphs[start_idx + 1:]:
        delete_paragraph(paragraph)

    pending_eq = None
    for block_type, text in BLOCKS:
        if block_type == "h1":
            add_heading(document, text, 1)
        elif block_type == "h2":
            add_heading(document, text, 2)
        elif block_type == "p":
            add_body(document, text)
        elif block_type == "eq":
            pending_eq = text
        elif block_type == "eq_no":
            if pending_eq is None:
                raise ValueError("公式编号前缺少公式内容。")
            add_equation(document, pending_eq, text)
            pending_eq = None
        elif block_type == "table":
            add_module_table(document, text)
        elif block_type == "figure":
            add_figure(document, text)
        else:
            raise ValueError(f"未知块类型: {block_type}")

    if pending_eq is not None:
        raise ValueError("存在未编号公式。")

    document.save(str(OUTPUT))

    section_text = []
    capture = False
    for paragraph in Document(str(OUTPUT)).paragraphs:
        text = paragraph.text.strip()
        if text == "四、专业实践总结报告":
            capture = True
            continue
        if capture and text:
            section_text.append(text)
    char_count = sum(len(t.replace(" ", "")) for t in section_text)
    print(f"Output: {OUTPUT}")
    print(f"Section-4 character count: {char_count}")


if __name__ == "__main__":
    main()
