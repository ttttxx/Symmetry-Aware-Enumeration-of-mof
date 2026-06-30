import os
import sys
import pytest

# 将项目根目录加入系统路径，以便能够导入你的模块
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from main import enrich_graph_data_with_topology, load_graph
from symmetry import SymmetryEngine
from enumerator import MOFEnumerator
from assembler import StructureAssembler
from verification import BridgeVerifier

# ================= 测试配置 =================
# 基于你提供的需求指定的固定参数
DATA_PREFIX = "ffc"
ATOM_CHOICES = ["4c_Cu_1_Ch"]
FORCED_EDGE = "1B_4H_Ch"


# ============================================

@pytest.fixture
def vba_graph_data():
    """
    Pytest Fixture：预先加载和处理 VBA 拓扑数据，供后续所有测试用例复用。
    """
    # 动态获取项目根目录下的 data 路径
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    graph_path = os.path.join(base_dir, "data", DATA_PREFIX, f"{DATA_PREFIX}_graph.json")

    assert os.path.exists(graph_path), f"找不到拓扑图文件: {graph_path}"

    # 模拟 main.py 中的数据加载和拓扑丰富过程
    graph_data = load_graph(graph_path)
    graph_data = enrich_graph_data_with_topology(graph_data)

    # 确保关键字段存在
    assert "alpha_matrix" in graph_data, "拓扑图中缺少 alpha_matrix，可能建图过程失败"
    assert "Bstar" in graph_data, "拓扑图中缺少 Bstar 矩阵"

    return graph_data


def test_symmetry_engine(vba_graph_data):
    """
    单元测试：测试 SymmetryEngine 能否正确找出自同构群
    """
    engine = SymmetryEngine(vba_graph_data, verbose=False)
    perms = engine.compute_automorphism_group()

    # VBA 拓扑应该具有对称性，即群大小应大于 0
    assert len(perms) > 0, "SymmetryEngine 未能找到任何对称性排列"
    # 确保恒等变换存在
    identity = tuple(range(vba_graph_data["n_nodes"]))
    assert identity in perms, "自同构群中缺失恒等变换 (Identity)"


def test_mof_enumerator(vba_graph_data):
    """
    单元测试：测试 MOFEnumerator 是否能根据节点和原子类型进行结构枚举
    """
    # 依赖对称性引擎的输出
    engine = SymmetryEngine(vba_graph_data, verbose=False)
    perms = engine.compute_automorphism_group()

    num_nodes = vba_graph_data["n_nodes"]
    enumerator = MOFEnumerator(perms, num_nodes, verbose=False)

    results = enumerator.enumerate(ATOM_CHOICES)

    assert len(results) > 0, "Enumerator 未能生成任何结构"

    # 验证输出的数据结构是否符合预期
    first_result = results[0]
    assert "coloring" in first_result, "枚举结果缺失 'coloring' 字段"
    assert "stabilizer_size" in first_result, "枚举结果缺失 'stabilizer_size' 字段"
    assert len(first_result["coloring"]) == num_nodes, "着色方案的节点数量与拓扑节点数不匹配"


def test_verification_logic(vba_graph_data):
    """
    单元测试：验证 BridgeVerifier 的 Burnside Lemma 自验证逻辑
    """
    engine = SymmetryEngine(vba_graph_data, verbose=False)
    perms = engine.compute_automorphism_group()

    # 只传入一个选项时，计算量小，理论可以完全穷举
    enumerator = MOFEnumerator(perms, vba_graph_data["n_nodes"], verbose=False)
    results = enumerator.enumerate(ATOM_CHOICES)

    verifier = BridgeVerifier(len(perms))
    is_verified, stats = verifier.verify(results, vba_graph_data["n_nodes"], ATOM_CHOICES, None)

    assert is_verified is True, "算法自验证 (Burnside Lemma) 失败"


def test_assembler_integration(vba_graph_data, tmp_path):
    """
    集成测试：结合组装器 (Assembler) 生成最终的 CIF 文件。
    使用 pytest 的 tmp_path 避免在测试时向真实的 results 文件夹写入垃圾文件。
    ffc 拓扑包含度-3 和度-4 节点，需要按度匹配配体。
    """
    # 1. 计算各拓扑节点的度，按度分配匹配的配体
    node_degrees = {}
    for e in vba_graph_data['edges']:
        s, t = str(e['source']), str(e['target'])
        node_degrees[s] = node_degrees.get(s, 0) + 1
        node_degrees[t] = node_degrees.get(t, 0) + 1

    structure_map = {}
    for i, node_info in enumerate(vba_graph_data['nodes']):
        node_id = str(node_info.get('id', i))
        deg = node_degrees.get(node_id, 0)
        if deg == 3:
            structure_map[node_id] = "3c_B_Ch"      # 3 个 X 连接点
        else:
            structure_map[node_id] = "4c_Cu_1_Ch"     # 4 个 X 连接点

    # 2. 初始化组装器并输出到临时目录
    assembler = StructureAssembler(output_dir=str(tmp_path))

    out_cif_filename = os.path.join(tmp_path, "test_ffc_assembled.cif")

    # 执行组装
    try:
        assembler.assemble_and_save(
            filename=out_cif_filename,
            structure_map=structure_map,
            graph_data=vba_graph_data,
            forced_edge=FORCED_EDGE
        )
    except Exception as e:
        pytest.fail(f"组装过程发生异常: {e}")

    # 3. 断言验证
    assert os.path.exists(out_cif_filename), "Assembler 未能生成 CIF 文件"
    assert os.path.getsize(out_cif_filename) > 0, "生成的 CIF 文件为空白"