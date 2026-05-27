# core/enumerator.py
import time
import itertools
import heapq
from collections import defaultdict
from typing import List, Tuple, Dict, Optional, Set


# ================= 辅助函数 (保持不变) =================
def cycles_to_perm(cycles: List[Tuple[int, ...]], n: int) -> List[int]:
    perm = list(range(n))
    for cycle in cycles:
        l = len(cycle)
        for k in range(l):
            perm[cycle[k]] = cycle[(k + 1) % l]
    return perm


def perm_to_cycles(perm: Tuple[int, ...], n: int) -> List[List[int]]:
    visited = [False] * n
    cycles = []
    for i in range(n):
        if not visited[i]:
            curr = i
            cycle = []
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = perm[curr]
            cycles.append(cycle)
    return cycles


# ================= Trie 节点定义 (保持不变) =================
class _TrieNode:
    __slots__ = ("children", "class_id", "stabilizer_size")

    def __init__(self):
        self.children: Dict[str, "_TrieNode"] = {}
        self.class_id: Optional[int] = None
        self.stabilizer_size: int = 0


# ================= 核心枚举器 =================
class MOFEnumerator:
    def __init__(self, permutations: List[Tuple[int, ...]], num_nodes: int,
                 verbose: bool = True,
                 topology_requirements: Optional[Dict[int, int]] = None,
                 user_selected_ligands: Optional[Dict[int, List[str]]] = None,
                 node_degrees: Optional[List[int]] = None):
        self.n = num_nodes
        self.perms = permutations
        self.perms_inv = self._compute_inverses(permutations)
        self.verbose = verbose
        self.topology_requirements = topology_requirements or {}
        self.user_selected_ligands = user_selected_ligands or {}
        self.node_degrees = node_degrees or []

        self._trie_root = _TrieNode()
        self.classes: Dict[int, Dict] = {}
        self.top_k_heap = []

    def _compute_inverses(self, perms):
        invs = []
        for p in perms:
            inv = [0] * self.n
            for i, j in enumerate(p):
                inv[j] = i
            invs.append(tuple(inv))
        return invs

    # ---------------- 基础工具 ----------------
    def _apply_perm(self, coloring: Tuple[str, ...], perm: Tuple[int, ...]) -> Tuple[str, ...]:
        n = self.n
        new_c = [""] * n
        for i in range(n):
            new_c[perm[i]] = coloring[i]
        return tuple(new_c)

    def _get_canonical(self, coloring: Tuple[str, ...]) -> Tuple[str, ...]:
        """获取正则形式"""
        canonical = coloring
        for p in self.perms:
            candidate = self._apply_perm(coloring, p)
            if candidate < canonical:
                canonical = candidate
        return canonical

    def _compute_stabilizer_and_orbit(self, coloring: Tuple[str, ...]) -> Tuple[int, List[Tuple[str, ...]]]:
        orbit_set = set()
        stab_count = 0
        for p in self.perms:
            mapped = self._apply_perm(coloring, p)
            orbit_set.add(mapped)
            if mapped == coloring:
                stab_count += 1
        return stab_count, list(orbit_set)

    def _compute_orbits_by_gen(self, generators: List[Tuple[int, ...]]) -> List[List[int]]:
        parent = list(range(self.n))

        def find(i):
            if parent[i] == i: return i
            parent[i] = find(parent[i])
            return parent[i]

        def union(i, j):
            root_i = find(i)
            root_j = find(j)
            if root_i != root_j: parent[root_i] = root_j

        for p in generators:
            for i, target in enumerate(p):
                union(i, target)

        orbit_dict = defaultdict(list)
        for i in range(self.n):
            orbit_dict[find(i)].append(i)
        return sorted(list(orbit_dict.values()), key=lambda x: x[0])

    # ---------------- 结果管理 ----------------
    def _trie_lookup(self, canonical: Tuple[str, ...]) -> Optional[int]:
        node = self._trie_root
        for c in canonical:
            if c not in node.children: return None
            node = node.children[c]
        return node.class_id

    def _trie_insert(self, canonical: Tuple[str, ...], class_id: int, stabilizer_size: int):
        node = self._trie_root
        for c in canonical:
            if c not in node.children: node.children[c] = _TrieNode()
            node = node.children[c]
        node.class_id = class_id
        node.stabilizer_size = stabilizer_size

    def _try_add_result(self, coloring: List[str], max_results: Optional[int]):
        raw_tuple = tuple(coloring)
        canonical_tup = self._get_canonical(raw_tuple)

        if self._trie_lookup(canonical_tup) is not None:
            return

        stab_size, _ = self._compute_stabilizer_and_orbit(canonical_tup)

        cid = len(self.classes) + 1000000
        should_store = True

        if max_results is not None:
            if len(self.classes) < max_results:
                heapq.heappush(self.top_k_heap, (stab_size, cid))
            else:
                min_stab, min_cid = self.top_k_heap[0]
                if stab_size > min_stab:
                    heapq.heapreplace(self.top_k_heap, (stab_size, cid))
                    if min_cid in self.classes: del self.classes[min_cid]
                else:
                    should_store = False

        self._trie_insert(canonical_tup, cid, stab_size)
        if should_store:
            self.classes[cid] = {
                "id": cid,
                "coloring": list(canonical_tup),
                "stabilizer_size": stab_size,
                "orbit_size": len(self.perms) // stab_size
            }

    # ---------------- [核心] 统一的 Lattice DFS (含早期剪枝) ----------------
    def _dfs_lattice_pruning(self, orbit_idx, orbits, atom_types, current_coloring, active_perm_indices, max_results,
                             state, current_degree_counts=None):
        state['steps'] += 1
        if state['steps'] > state['limit']:
            return

        # 初始化按度数计数的计数器 (首次调用)
        if current_degree_counts is None and self.topology_requirements:
            current_degree_counts = {deg: 0 for deg in self.topology_requirements}

        # Base Case
        if orbit_idx == len(orbits):
            self._try_add_result(current_coloring, max_results)
            return

        current_nodes = orbits[orbit_idx]
        # 获取该轨道中节点的度数 (同一轨道所有节点度数相同)
        orbit_degree = self.node_degrees[current_nodes[0]] if self.node_degrees else None

        # 获取该度数可用的配体列表 — 从 user_selected_ligands 动态选取
        tryable_atoms = atom_types
        if orbit_degree is not None and orbit_degree in self.user_selected_ligands:
            tryable_atoms = self.user_selected_ligands[orbit_degree]

        for atom in tryable_atoms:
            for node in current_nodes:
                current_coloring[node] = atom

            # ---- 剪枝 1: 按度数统计的总量约束检查 ----
            if self.topology_requirements and current_degree_counts is not None and orbit_degree is not None:
                new_count = current_degree_counts[orbit_degree] + len(current_nodes)
                if new_count > self.topology_requirements.get(orbit_degree, new_count):
                    for node in current_nodes:
                        current_coloring[node] = None
                    continue  # 该度数节点已超标，剪掉此分支
                current_degree_counts[orbit_degree] = new_count

            # ---- 剪枝 2: 正则对称性剪枝 ----
            next_active_indices = []
            is_pruned = False

            for p_idx in active_perm_indices:
                p_inv = self.perms_inv[p_idx]
                g_is_smaller = False
                g_is_larger = False

                for i in range(self.n):
                    val_i = current_coloring[i]
                    if val_i is None: break

                    src = p_inv[i]
                    val_src = current_coloring[src]
                    if val_src is None: break

                    if val_src < val_i:
                        g_is_smaller = True
                        break
                    elif val_src > val_i:
                        g_is_larger = True
                        break

                if g_is_smaller:
                    is_pruned = True
                    break
                if not g_is_larger:
                    next_active_indices.append(p_idx)

            if not is_pruned:
                self._dfs_lattice_pruning(orbit_idx + 1, orbits, atom_types, current_coloring,
                                          next_active_indices, max_results, state, current_degree_counts)
                if state['steps'] > state['limit']:  # 熔断返回
                    for node in current_nodes: current_coloring[node] = None
                    return

            # Backtrack 按度数计数器
            if self.topology_requirements and current_degree_counts is not None and orbit_degree is not None:
                current_degree_counts[orbit_degree] -= len(current_nodes)

            # Backtrack
            for node in current_nodes:
                current_coloring[node] = None

    # ---------------- 策略：Unified Lattice Descent ----------------
    def enumerate_subgroups(self, atom_types: List[str], max_results: int = None):
        # 任务调度堆（最小堆，按 score 升序弹出）
        # 注意：与结果追踪堆 self.top_k_heap 相互独立，互不干扰
        task_heap = []
        _task_idx = 0  # tie-breaker，避免字典直接比较引发 TypeError

        # 1. 全群 (G)
        orbits_g = self._compute_orbits_by_gen(self.perms)
        heapq.heappush(task_heap, (len(orbits_g), _task_idx, {"name": "Full Group (G)", "orbits": orbits_g, "score": len(orbits_g)}))
        _task_idx += 1

        # 2. 循环子群 (Cyclic)
        seen_fp = set()
        for p in self.perms:
            orbits = self._compute_orbits_by_gen([p])
            fp = tuple(sorted([len(o) for o in orbits]))
            if fp not in seen_fp:
                seen_fp.add(fp)
                heapq.heappush(task_heap, (len(orbits), _task_idx, {"name": f"Cyclic Type {fp[:2]}...", "orbits": orbits, "score": len(orbits)}))
                _task_idx += 1

        # 3. 恒等群 (Identity) - 兜底
        identity_orbits = [[i] for i in range(self.n)]
        heapq.heappush(task_heap, (self.n, _task_idx, {"name": "Identity (Fallback)", "orbits": identity_orbits, "score": self.n}))
        _task_idx += 1

        if self.verbose:
            print(f"  [Lattice] Generated {len(task_heap)} tasks. Executing Best-N Strategy...")

        # 4. 执行
        initial_active = list(range(1, len(self.perms)))

        # [修改] 重新引入严格的 Space Check 和 适度的 Step Limit
        SPACE_LIMIT = 5 * 10 ** 6  # 500万：超过这个大小的中间任务直接跳过
        DFS_STEP_LIMIT = 200000  # 20万：单次 DFS 最多跑这么多步 (约 1-2秒)

        while task_heap:
            _, _, task = heapq.heappop(task_heap)  # 按 score 升序弹出（高对称优先）
            # 如果已经找够了，且当前不是必须跑的高对称任务，可以提前观察退出
            if max_results and len(self.classes) >= max_results:
                # 只有当任务非常难算时才退出，否则为了 Best-N 我们可以多跑跑
                pass

            name = task['name']
            orbits = task['orbits']
            n_orbits = len(orbits)
            space_size = len(atom_types) ** n_orbits

            is_identity = (name == "Identity (Fallback)")

            # [熔断器逻辑 — 对应论文 Algorithm 1 第10行]
            # Tier 1 (全群) / Tier 2 (循环子群)：不受 SPACE_LIMIT 约束，仅受 DFS_STEP_LIMIT 步数保护
            # Tier 3 (恒等群/Identity)：若预测状态空间 > SPACE_LIMIT 则熔断跳过
            if is_identity and space_size > SPACE_LIMIT:
                if self.verbose:
                    print(f"    [Skip] {name} (Tier 3, Space {space_size:.1e} > {SPACE_LIMIT:.1e})")
                continue

            # 打印日志防止用户以为卡死
            if self.verbose:
                print(f"    [Run] {name} (Orbits: {n_orbits})...")

            current_coloring = [None] * self.n
            state = {'steps': 0, 'limit': DFS_STEP_LIMIT}

            self._dfs_lattice_pruning(0, orbits, atom_types, current_coloring, initial_active, max_results, state)

            # 如果是 Identity 并且跑满了步数，说明真的很难找了
            if is_identity and state['steps'] > state['limit']:
                if self.verbose: print(f"    [Stop] {name} hit step limit.")

    # ---------------- 主入口 ----------------
    def enumerate(self, atom_types: List[str], composition: Dict[str, int] = None, max_results: int = None):
        self._trie_root = _TrieNode()
        self.classes = {}
        self.top_k_heap = []

        # 统一执行
        self.enumerate_subgroups(atom_types, max_results)

        # 排序输出
        res = sorted(list(self.classes.values()), key=lambda x: x['stabilizer_size'], reverse=True)
        for i, r in enumerate(res): r['id'] = i

        self._verify_burnside()
        return res

    def _verify_burnside(self):
        lhs = sum(1.0 / d['stabilizer_size'] for d in self.classes.values())
        print(f"  [Verification] Sum(1/|Stab|) = {lhs:.4f}")