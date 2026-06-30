"""
CIF 文件名解析模块。

解析规则：
  - 去掉 .cif 后缀，按 _v\\d+- 分割
  - 第一段 = 拓扑名
  - 后续以 vN- 开头的段 = 节点列表（数量不固定）
  - 最后一段中，去除已知 edge 后缀后剩余部分 = 最后一个节点名
  - edge 名通过 edges/ 目录动态发现
"""

import os
import re

_EDGES_CACHE = None


def _load_known_edges(edges_dir=None):
    """从 edges/ 目录加载所有已知 edge 名（按长度降序，确保最长优先匹配）。"""
    global _EDGES_CACHE
    if _EDGES_CACHE is not None:
        return _EDGES_CACHE
    if edges_dir is None:
        edges_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "edges")
    edges = []
    if os.path.isdir(edges_dir):
        for f in os.listdir(edges_dir):
            if f.endswith(".cif"):
                edges.append(f[:-4])
    edges.sort(key=len, reverse=True)
    _EDGES_CACHE = edges
    return _EDGES_CACHE


def parse_filename(filename, known_edges=None):
    """
    解析 CIF 文件名，返回 (topology: str, nodes: list[str], edge: str)。

    >>> parse_filename("cdl_v1-4c_Cu_1_Ch_1B_4H_Ch.cif", ["1B_4H_Ch"])
    ('cdl', ['4c_Cu_1_Ch'], '1B_4H_Ch')

    >>> parse_filename("ffc_v1-4c_Cu_1_Ch_v2-4c_Cu_1_Ch_v3-3c_B_Ch_1B_4H_Ch.cif", ["1B_4H_Ch"])
    ('ffc', ['4c_Cu_1_Ch', '4c_Cu_1_Ch', '3c_B_Ch'], '1B_4H_Ch')
    """
    if known_edges is None:
        known_edges = _load_known_edges()

    name = filename
    if name.endswith(".cif"):
        name = name[:-4]

    # 按 _v\\d+- 分割，保留节点边界
    parts = re.split(r"_(?=v\d+-)", name)

    topology = parts[0]

    # 提取每个 vN- 段的内容（去掉 vN- 前缀）
    raw_nodes = []
    for part in parts[1:]:
        m = re.match(r"v\d+-(.+)", part)
        if m:
            raw_nodes.append(m.group(1))

    if not raw_nodes:
        raise ValueError(f"未找到任何节点标记 (vN-): {filename}")

    # 最后一段可能含 edge 后缀，尝试匹配拆分
    edge = ""
    last_node = raw_nodes[-1]
    for known_edge in known_edges:
        suffix = "_" + known_edge
        if last_node.endswith(suffix):
            raw_nodes[-1] = last_node[: -len(suffix)]
            edge = known_edge
            break

    if not edge:
        raise ValueError(
            f"无法从最后一段中识别 edge 名。最后一段='{last_node}'，"
            f"已知 edge={known_edges[:5]}..."
        )

    return topology, raw_nodes, edge


# ── 单元测试 ──

def _test():
    test_edges = ["1B_4H_Ch", "1B_1TrU"]

    # 1 节点
    topo, nodes, edge = parse_filename("cdl_v1-4c_Cu_1_Ch_1B_4H_Ch.cif", test_edges)
    assert topo == "cdl", f"topo: {topo}"
    assert nodes == ["4c_Cu_1_Ch"], f"nodes: {nodes}"
    assert edge == "1B_4H_Ch", f"edge: {edge}"

    # 2 节点
    topo, nodes, edge = parse_filename("fog_v1-3c_B_Ch_v2-4c_Cu_1_Ch_1B_4H_Ch.cif", test_edges)
    assert topo == "fog"
    assert nodes == ["3c_B_Ch", "4c_Cu_1_Ch"]
    assert edge == "1B_4H_Ch"

    # 3 节点
    topo, nodes, edge = parse_filename(
        "ffc_v1-4c_Cu_1_Ch_v2-4c_Cu_1_Ch_v3-3c_B_Ch_1B_4H_Ch.cif", test_edges
    )
    assert topo == "ffc"
    assert nodes == ["4c_Cu_1_Ch", "4c_Cu_1_Ch", "3c_B_Ch"]
    assert edge == "1B_4H_Ch"

    # 3 节点含 triazine
    topo, nodes, edge = parse_filename(
        "fjh_v1-3c_triazine_Ch_v2-4c_Cu_1_Ch_v3-4c_Cu_1_Ch_1B_4H_Ch.cif",
        test_edges,
    )
    assert topo == "fjh"
    assert nodes == ["3c_triazine_Ch", "4c_Cu_1_Ch", "4c_Cu_1_Ch"]
    assert edge == "1B_4H_Ch"

    # 4 节点（合成用例）
    topo, nodes, edge = parse_filename(
        "abcd_v1-4c_Cu_1_Ch_v2-3c_B_Ch_v3-4c_Cu_1_Ch_v4-3c_triazine_Ch_1B_4H_Ch.cif",
        test_edges,
    )
    assert topo == "abcd"
    assert nodes == ["4c_Cu_1_Ch", "3c_B_Ch", "4c_Cu_1_Ch", "3c_triazine_Ch"]
    assert edge == "1B_4H_Ch"

    # 使用真实 edges 目录的 end-to-end 测试
    topo, nodes, edge = parse_filename("nbo_v1-4c_Cu_1_Ch_1B_4H_Ch.cif")
    assert topo == "nbo"
    assert nodes == ["4c_Cu_1_Ch"]
    assert edge == "1B_4H_Ch"

    print("全部 6 项测试通过 ✓")


if __name__ == "__main__":
    _test()