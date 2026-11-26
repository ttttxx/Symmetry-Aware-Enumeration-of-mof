import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const TMP = "/private/tmp/codex-presentations/20260712-mof-final";
const SOURCE = "/Users/jing/Desktop/专业实践（中期）/专业实践中期答辩-田婧.pptx";
const STARTER = path.join(TMP, "template-starter.pptx");
const FINAL = path.join(TMP, "专业实践终期答辩-田婧.pptx");
const MEDIA = path.join(TMP, "manual-media");
const MAP = path.join(TMP, "template-frame-map.json");

const edit = (ids, action = "rewrite") => ids.map((sourceElementId) => ({ sourceElementId, action }));

const outputSlides = [
  { outputSlide: 1, sourceSlide: 1, narrativeRole: "opening thesis", reuseMode: "duplicate-slide", editTargets: edit(["sh/ts7md4r2", "sh/sryl4zqx", "sh/ozmd8r65", "sh/mx4v6hoz", "sh/crulkr2h"]) },
  { outputSlide: 2, sourceSlide: 2, narrativeRole: "agenda", reuseMode: "duplicate-slide", editTargets: edit(["sh/ru5cbqd4", "sh/a10jqpsj", "sh/bip8jmho", "sh/65kfmpor"]) },
  { outputSlide: 3, sourceSlide: 3, narrativeRole: "background divider", reuseMode: "duplicate-slide", editTargets: [] },
  { outputSlide: 4, sourceSlide: 4, narrativeRole: "background and stakes", reuseMode: "duplicate-slide", editTargets: edit(["sh/5cva1cfq", "sh/sna103ap"]) },
  { outputSlide: 5, sourceSlide: 5, narrativeRole: "progress divider", reuseMode: "duplicate-slide", editTargets: [] },
  { outputSlide: 6, sourceSlide: 6, narrativeRole: "completed plan", reuseMode: "duplicate-slide", editTargets: edit(["sh/zi98nu94", "sh/h4bupgn6", "sh/gju58ji9", "sh/1k3mhozu"]) },
  { outputSlide: 7, sourceSlide: 7, narrativeRole: "three-stage progression", reuseMode: "duplicate-slide", editTargets: edit(["sh/dcbud0ra", "sh/e9c7ixkv", "sh/eh4jil4r", "sh/id4ju1oz", "sh/ozq5szq9", "sh/q143ylov", "sh/f6xkvq5c"]) },
  { outputSlide: 8, sourceSlide: 15, narrativeRole: "final-results divider", reuseMode: "duplicate-slide", editTargets: edit(["sh/snq9sba5"]) },
  { outputSlide: 9, sourceSlide: 16, narrativeRole: "supporting engineering results", reuseMode: "duplicate-slide", editTargets: edit(["sh/e1sf65o3", "sh/fm1gzq5o", "sh/1ojy10ne", "sh/ih0jedob", "sh/tcbmdcre", "sh/ji9k7ypw", "sh/sb2l47at"]) },
  { outputSlide: 10, sourceSlide: 20, narrativeRole: "MOF system thesis", reuseMode: "duplicate-slide", editTargets: edit(["sh/po3itc72", "sh/qh4nupg7", "sh/grm5ofy5", "sh/rid4nuhs", "sh/32lgjid0", "sh/ulwz6d4b", "sh/gnyh8nm1"]) },
  { outputSlide: 11, sourceSlide: 22, narrativeRole: "system screenshot chain", reuseMode: "duplicate-slide", editTargets: [
    ...edit(["sh/gf65w7ql", "sh/qtwz2547", "sh/t8fah0f6", "sh/dkzi1gz2", "sh/cjqhsbih", "sh/bih0z6hc", "sh/qh8zq10r", "sh/xg3id0zm", "sh/fadobits", "sh/gbmpk3ud"]),
    ...edit(["im/ozex8ni9", "im/fa5gbihc", "im/e9cf2d07", "im/1cnyd8zi"], "replace"),
  ] },
  { outputSlide: 12, sourceSlide: 23, narrativeRole: "system overview", reuseMode: "duplicate-slide", editTargets: [
    ...edit(["sh/xovqhcz2", "sh/ex0filc3", "sh/hg7uhcra", "sh/xs7upszy", "sh/ju9wr2ho", "sh/6h0v2x0f", "sh/5obulcfy", "sh/ux0bmxwb", "sh/8vitknel", "sh/fe5cf294", "sh/p47elcr6", "sh/6l8neh8b", "sh/snq5grq1", "sh/hsj6dw7y", "sh/vux4z29o"]),
    ...edit(["im/9gjalwra"], "replace"),
  ] },
  { outputSlide: 13, sourceSlide: 25, narrativeRole: "manual evidence gallery", reuseMode: "duplicate-slide", editTargets: [
    ...edit(["sh/gbq1obep", "sh/5kfi9cf2", "sh/4j61g7eh", "sh/hcvylsb6", "sh/mlkfmdcj", "sh/9obyh8vu", "sh/yh0fitc7", "sh/q9kve9g3", "sh/dcbup4fe", "sh/250bapgr", "sh/p8bul4fi", "sh/uh0b6pgv", "sh/kny9gvet"]),
    ...edit(["im/7md03m94", "im/6l4juh8j", "im/2h832lcr", "im/3ihkvqtc", "im/4jql4bah", "im/pkz2xgb2"], "replace"),
  ] },
  { outputSlide: 14, sourceSlide: 27, narrativeRole: "end-to-end workflow", reuseMode: "duplicate-slide", editTargets: [
    ...edit(["sh/utobax4f", "sh/ulgzuh8r", "sh/rit0vqtk", "sh/aho3y9cj", "sh/l8zmxov6", "sh/5836tsfa", "sh/i5c7yxwz", "sh/u1s7udwn", "sh/bah8byh4", "sh/xwzqd8zu", "sh/jy18fih0"]),
    { sourceElementId: "sh/kn694nmx", action: "delete" },
  ] },
  { outputSlide: 15, sourceSlide: 21, narrativeRole: "implementation difficulties", reuseMode: "duplicate-slide", editTargets: edit(["sh/ra9kryp8", "sh/psfator6", "sh/1kru5kj2", "sh/2l0vypkn", "sh/y9cfatgv", "sh/zalgjyhg", "sh/14vupk7a", "sh/yxov69on", "sh/10ba50na", "sh/gf29cvmp"]) },
  { outputSlide: 16, sourceSlide: 17, narrativeRole: "difficulty divider", reuseMode: "duplicate-slide", editTargets: edit(["sh/gr6hkv6l"]) },
  { outputSlide: 17, sourceSlide: 18, narrativeRole: "next steps", reuseMode: "duplicate-slide", editTargets: edit(["sh/hsn2l4bu", "sh/0netwrex", "sh/v2xovu9o", "sh/je5knut0"]) },
  { outputSlide: 18, sourceSlide: 19, narrativeRole: "closing", reuseMode: "duplicate-slide", editTargets: edit(["sh/xkbq10be", "sh/wj2psvat", "sh/wfy1s7ix", "sh/yhcfydgj"]) },
];

const frameMap = {
  outputSlides,
  omittedSourceSlides: [
    { sourceSlide: 8, reason: "detailed virtual-machine template pages omitted; the supporting engineering summary is retained" },
    { sourceSlide: 9, reason: "detailed virtual-machine template pages omitted; the supporting engineering summary is retained" },
    { sourceSlide: 10, reason: "detailed MooseFS pages omitted; the supporting engineering summary is retained" },
    { sourceSlide: 11, reason: "detailed MooseFS pages omitted; the supporting engineering summary is retained" },
    { sourceSlide: 12, reason: "detailed testing-tool pages omitted; the supporting engineering summary is retained" },
    { sourceSlide: 13, reason: "detailed testing-tool pages omitted; the supporting engineering summary is retained" },
    { sourceSlide: 14, reason: "detailed testing-tool pages omitted; the supporting engineering summary is retained" },
    { sourceSlide: 24, reason: "unused six-image layout; replaced by the selected manual screenshot gallery" },
    { sourceSlide: 26, reason: "duplicate three-stage layout; the final deck uses the earlier workflow slide" },
    { sourceSlide: 28, reason: "unused four-stage layout; final deck closes with the source closing slide" },
  ],
};

const bytes = async (file) => {
  const b = await fs.readFile(file);
  return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);
};

const imageContentType = (file) => {
  const ext = path.extname(file).toLowerCase();
  return ext === ".jpg" || ext === ".jpeg" ? "image/jpeg" : ext === ".svg" ? "image/svg+xml" : "image/png";
};

function shapeById(slide, id) {
  const shape = slide.shapes.items.find((item) => String(item.id) === String(id));
  if (!shape) throw new Error(`Missing shape ${id} on slide ${slide.slideNumber ?? "?"}`);
  return shape;
}

function imageById(slide, id) {
  const image = slide.images.items.find((item) => String(item.id) === String(id));
  if (!image) throw new Error(`Missing image ${id} on slide ${slide.slideNumber ?? "?"}`);
  return image;
}

function setText(slide, id, text) {
  shapeById(slide, id).text = text;
}

function deleteShape(slide, id) {
  const shape = shapeById(slide, id);
  if (typeof shape.delete === "function") shape.delete();
  else shape.text = "";
}

async function replaceImage(slide, id, filename, frame, alt) {
  const image = imageById(slide, id);
  await image.replace({ blob: await bytes(path.join(MEDIA, filename)), contentType: imageContentType(filename), alt });
  image.frame = frame;
  image.fit = "contain";
}

async function writeMap() {
  await fs.writeFile(MAP, `${JSON.stringify(frameMap, null, 2)}\n`, "utf8");
}

async function editDeck() {
  const presentation = await PresentationFile.importPptx(await FileBlob.load(STARTER));

  // Source slide 1: update identity fields while retaining the title-slide composition.
  let s = presentation.slides.getItem(0);
  setText(s, 6, "专业实践终期答辩");
  setText(s, 7, "Professional Practice Final Defense");
  setText(s, 20, "实践基地");
  setText(s, 25, "实践基地");
  setText(s, 27, "航天软件评测中心 / 北京化工大学");
  setText(s, 223, "2026 / 07 / 12");

  // Agenda and front-half framing.
  s = presentation.slides.getItem(1);
  setText(s, 40, "项目背景与目标");
  setText(s, 5, "前期工程实践");
  setText(s, 17, "MOF系统终期成果");
  setText(s, 23, "难点与后续计划");

  s = presentation.slides.getItem(3);
  setText(s, 11, "项目背景与目标");
  setText(s, 7, "本项目由前期工程环境建设和后期科研软件开发两部分组成。前期围绕国产化环境、MooseFS 分布式存储和测试工具链开展实践；后期将 MOF 对称性感知结构枚举流程转化为可操作的浏览器端系统，形成从拓扑配置到结构导出的完整链路。终期工作的核心目标，是让算法结果能够被研究人员直接配置、查看、复核和归档。");

  s = presentation.slides.getItem(5);
  setText(s, 3, "实践计划的完成情况");
  setText(s, 10, "工程基础与存储");
  setText(s, 28, "MOF 枚举与验证");
  setText(s, 29, "3D 展示与导出");

  s = presentation.slides.getItem(6);
  setText(s, 3, "三段递进式终期工作");
  setText(s, 19, "工程基础");
  setText(s, 24, "国产化环境、MooseFS 集群与测试工具链，建立可运行的实践基础");
  setText(s, 31, "核心系统");
  setText(s, 29, "完成拓扑需求分析、对称性感知枚举、结果验证与 CIF 结构生成");
  setText(s, 39, "终期交付");
  setText(s, 23, "完成 Streamlit 界面、py3Dmol 三维查看、结构汇总、CIF/JSON 导出及配套说明书");
  setText(s, 73, "阶段一：2025/3—2025/8\n阶段二：2025/9—2026/3\n阶段三：2026/4—2026/6");
  setText(s, 20, "阶段一：2025/3—2025/8\n阶段二：2025/9—2026/3\n阶段三：2026/4—2026/6");

  // Section 3 divider and supporting results.
  s = presentation.slides.getItem(7);
  setText(s, 87, "终期成果与系统价值");

  s = presentation.slides.getItem(8);
  setText(s, 11, "前期工程实践为系统开发提供支撑");
  setText(s, 10, "环境标准化");
  setText(s, 12, "完成银河麒麟 V10 模板与软件仓库管理，为后续开发提供可复用环境。");
  setText(s, 28, "集群高可用");
  setText(s, 30, "参与 MooseFS 集群高可用配置与测试，积累数据可靠性和故障处理经验。");
  setText(s, 29, "工具链优化");
  setText(s, 31, "完成日志分析和嵌入式测试度量功能，为终期软件交付建立工程基础。");

  // MOF system thesis slide.
  s = presentation.slides.getItem(9);
  setText(s, 5, "MOF 系统聚焦两个问题");
  setText(s, 12, "算法如何落地");
  setText(s, 19, "Part 01");
  setText(s, 13, "将拓扑、配体、枚举、验证和 CIF 生成组织成可操作的 Web 流程，降低命令行与外部软件之间的切换成本。");
  setText(s, 42, "结果如何复核");
  setText(s, 40, "Part 02");
  setText(s, 38, "用对称性指标、结构表格、3D 视图以及 CIF/JSON 导出建立可复核链路；系统用于候选筛选，不直接替代稳定性计算。");

  // Manual screenshot chain.
  s = presentation.slides.getItem(10);
  setText(s, 15, "系统功能链：配置—计算—审查—归档");
  setText(s, 31, "说明书中的真实界面截图，展示一次运行从参数输入到结果导出的关键节点");
  setText(s, 90, "界面配置");
  setText(s, 104, "启动系统并设置拓扑、配体和枚举参数");
  setText(s, 105, "动态需求");
  setText(s, 106, "根据拓扑节点度数自动生成配体选项");
  setText(s, 107, "三维审查");
  setText(s, 108, "并行查看 Top-K 候选结构和渲染样式");
  setText(s, 111, "数据归档");
  setText(s, 112, "逐结构下载 CIF，并导出完整 JSON");
  await replaceImage(s, 113, "image9.png", { left: 82.18, top: 285.27, width: 252.69, height: 176.15 }, "MOF 系统主界面");
  await replaceImage(s, 114, "image10.png", { left: 371.1, top: 285.27, width: 252.69, height: 176.15 }, "拓扑动态配体配置");
  await replaceImage(s, 115, "image12.png", { left: 660.02, top: 285.27, width: 252.69, height: 176.15 }, "三维多视口结构查看");
  await replaceImage(s, 116, "image13.png", { left: 948.95, top: 285.27, width: 252.69, height: 176.15 }, "CIF 与 JSON 数据导出");

  // Single-image system overview.
  s = presentation.slides.getItem(11);
  setText(s, 6, "系统主界面把参数、计算和结果放在同一条链路上");
  setText(s, 85, "系统主界面");
  setText(s, 91, "参数—结果");
  await replaceImage(s, 111, "image9.png", { left: 71.36, top: 118.38, width: 478.05, height: 528.95 }, "MOF 可视化系统主界面截图");
  setText(s, 42, "拓扑与配体");
  setText(s, 44, "动态读取拓扑需求，按节点度数提供配体选择。");
  setText(s, 49, "枚举与验证");
  setText(s, 51, "调用对称性分析和结构枚举流程，保留运行状态与日志。");
  setText(s, 56, "三维结果");
  setText(s, 58, "以指标卡、Top-K 视图和结构表同时呈现结果。");
  setText(s, 93, "可复核导出");
  setText(s, 95, "候选结构可逐个导出 CIF，完整结果可导出 JSON。");
  setText(s, 100, "真实使用场景");
  setText(s, 102, "研究人员可在浏览器中完成配置、审查、筛选和归档。");
  setText(s, 107, "工程边界");
  setText(s, 109, "结果用于候选结构准备，不替代独立稳定性计算。");

  // Six-image evidence gallery.
  s = presentation.slides.getItem(12);
  setText(s, 21, "软件说明书中的真实截图对应系统的六个关键环节");
  const gallery = [
    [28, 29, "主界面", "image3.jpeg", "启动后的参数配置区"],
    [37, 38, "动态配置", "image10.png", "拓扑切换后生成配体控件"],
    [41, 42, "运行状态", "image11.png", "计算过程与运行日志"],
    [45, 46, "三维审查", "image12.png", "Top-K 多视口渲染"],
    [49, 50, "CIF下载", "image6.png", "视口级结构文件导出"],
    [33, 34, "JSON归档", "image7.png", "完整结果预览与导出"],
  ];
  const galleryFrames = [
    { left: 192.48, top: 144.15, width: 181.42, height: 102.05 },
    { left: 546.25, top: 144.15, width: 181.42, height: 102.05 },
    { left: 900.02, top: 144.15, width: 181.42, height: 102.05 },
    { left: 192.48, top: 416.56, width: 181.42, height: 102.05 },
    { left: 546.25, top: 416.56, width: 181.42, height: 102.05 },
    { left: 900.02, top: 416.56, width: 181.42, height: 102.05 },
  ];
  const galleryImageIds = [58, 59, 60, 61, 62, 63];
  for (let i = 0; i < gallery.length; i++) {
    const [labelId, descId, label, file, desc] = gallery[i];
    setText(s, labelId, label);
    setText(s, descId, desc);
    await replaceImage(s, galleryImageIds[i], file, galleryFrames[i], desc);
  }

  // End-to-end flow.
  s = presentation.slides.getItem(13);
  setText(s, 8, "系统流程：拓扑输入—可复核输出");
  setText(s, 19, "配置输入");
  setText(s, 21, "选择拓扑、节点配体、边缘分子与 Top-K 参数");
  setText(s, 33, "结构生成");
  setText(s, 34, "分析对称性，枚举不等价着色并组装 CIF");
  setText(s, 41, "结果审查");
  setText(s, 42, "查看指标、3D 结构、汇总表和运行日志");
  setText(s, 46, "输出链路");
  setText(s, 50, "3D 视图：旋转、缩放、切换渲染样式");
  setText(s, 52, "结构表：按稳定子与轨道属性浏览");
  setText(s, 54, "文件：逐结构导出 CIF，汇总导出 JSON");
  deleteShape(s, 2);
  deleteShape(s, 29);

  // Difficulties and solutions.
  s = presentation.slides.getItem(14);
  setText(s, 4, "系统难点：多个环节必须同时成立");
  setText(s, 150, "算法正确只是起点：候选空间控制、几何组装、结果验证和交互式交付共同决定系统是否真正可用。");
  setText(s, 24, "枚举规模");
  setText(s, 23, "拓扑和配体组合会快速放大搜索空间，需要用对称性排序和上限控制保留有价值的候选。");
  setText(s, 41, "几何匹配");
  setText(s, 40, "抽象着色结果还要放回晶体坐标；节点、边缘分子与拓扑方向不匹配时会导致组装失败。");
  setText(s, 68, "结果可信");
  setText(s, 67, "结构数量不足以说明质量，需要结合稳定子、轨道和 Burnside 验证检查结果。");
  setText(s, 77, "交互落地");
  setText(s, 76, "计算、日志、3D、表格与文件导出需要稳定衔接。");

  // Difficulty divider and follow-up plan.
  s = presentation.slides.getItem(15);
  setText(s, 87, "难点与后续计划");
  s = presentation.slides.getItem(16);
  setText(s, 3, "难点对应的后续工作");
  setText(s, 40, "规模与稳定性\n继续验证大规模枚举、批量任务和复杂拓扑下的运行边界。");
  setText(s, 25, "可用性与验证\n完善几何质量检查、结果检索、异常提示和面向科研用户的交互体验。");
  setText(s, 5, "当前边界：系统面向候选结构生成、可视化和数据管理，后续仍需与独立计算流程衔接。");

  // Closing slide.
  s = presentation.slides.getItem(17);
  setText(s, 6, "感谢老师的悉心指导");
  setText(s, 7, "Professional Practice Final Defense");
  setText(s, 21, "THANKS");
  setText(s, 34, "2026 / 07 / 12");

  await fs.mkdir(path.dirname(FINAL), { recursive: true });
  const exported = await PresentationFile.exportPptx(presentation);
  await exported.save(FINAL);
  console.log(`Wrote ${FINAL}`);
}

await writeMap();
if (process.env.MAP_ONLY === "1") process.exit(0);
await editDeck();
