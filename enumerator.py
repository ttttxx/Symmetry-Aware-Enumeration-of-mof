# core/enumerator.py
import heapq
from collections import defaultdict
from typing import List, Tuple, Dict, Optional, Set


# ── Identity 子群枚举阈值配置 ──────────────────────────────────────────────
# 当估算不等价轨道数 (naive_space / |G|) 超过此值时，切换为限额枚举模式
IDENTITY_ORBIT_LIMIT: int = 500_000
# 限额模式下 Identity 任务最多新增的不等价结构数
IDENTITY_NEW_LIMIT: int = 10_000

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

    def _try_add_result(self, coloring: List[str]):
        raw_tuple = tuple(coloring)
        canonical_tup = self._get_canonical(raw_tuple)

        if self._trie_lookup(canonical_tup) is not None:
            return

        stab_size, _ = self._compute_stabilizer_and_orbit(canonical_tup)
        cid = len(self.classes)

        self._trie_insert(canonical_tup, cid, stab_size)
        self.classes[cid] = {
            "id": cid,
            "coloring": list(canonical_tup),
            "stabilizer_size": stab_size,
            "orbit_size": len(self.perms) // stab_size
        }

    # ---------------- [核心] 统一的 Lattice DFS (含早期剪枝) ----------------
    def _dfs_lattice_pruning(self, orbit_idx, orbits, atom_types, current_coloring, active_perm_indices,
                             current_degree_counts=None,
                             identity_start: int = None,
                             identity_budget: int = None):
        # 初始化按度数计数的计数器 (首次调用)
        if current_degree_counts is None and self.topology_requirements:
            current_degree_counts = {deg: 0 for deg in self.topology_requirements}

        # Base Case
        if orbit_idx == len(orbits):
            # 预算检查：限额模式下超额则跳过
            if identity_start is not None and identity_budget is not None:
                if len(self.classes) - identity_start >= identity_budget:
                    return
            self._try_add_result(current_coloring)
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
                                          next_active_indices, current_degree_counts,
                                          identity_start, identity_budget)

            # Backtrack 按度数计数器
            if self.topology_requirements and current_degree_counts is not None and orbit_degree is not None:
                current_degree_counts[orbit_degree] -= len(current_nodes)

            # Backtrack
            for node in current_nodes:
                current_coloring[node] = None

    # ---------------- 策略：Unified Lattice Descent ----------------
    def enumerate_subgroups(self, atom_types: List[str]):
        self._last_atom_types = atom_types  # 供 _verify_burnside 计算覆盖率
        # 任务调度堆（最小堆，按 orbit count 升序弹出，即最高对称性优先）
        task_heap = []
        _task_idx = 0  # tie-breaker，避免字典直接比较引发 TypeError

        # 1. 全群 (G) — 子群格顶端
        orbits_g = self._compute_orbits_by_gen(self.perms)
        heapq.heappush(task_heap, (len(orbits_g), _task_idx, {"name": "Full Group (G)", "orbits": orbits_g, "score": len(orbits_g)}))
        _task_idx += 1

        # 2. 完整子群格 (Lattice) — 中间真子群
        try:
            from sympy.combinatorics import Permutation, PermutationGroup

            def _enumerate_all_sympy_subgroups(G):
                """Bottom-up subgroup lattice enumeration via SymPy (no .subgroups() required)."""
                g_order = G.order()
                seen_keys = set()
                result = []

                def add_if_new(H):
                    key = frozenset(tuple(e.array_form) for e in H.generate())
                    if key not in seen_keys:
                        seen_keys.add(key)
                        result.append(H)
                        return True
                    return False

                # Seed with cyclic subgroups from every element
                for elem in G.generate():
                    add_if_new(PermutationGroup(elem))

                # Iteratively combine pairs until stable
                i = 0
                while i < len(result):
                    H1 = result[i]
                    for H2 in result[:i]:
                        K = PermutationGroup(*(list(H1.generators) + list(H2.generators)))
                        if K.order() < g_order:
                            add_if_new(K)
                    i += 1

                return result

            def _sympy_perm_to_tuple(p, n):
                arr = list(p.array_form)
                arr += list(range(len(arr), n))  # 补齐 SymPy 截断的末尾不动点
                return tuple(arr)

            sympy_perms = [Permutation(list(p), size=self.n) for p in self.perms]
            G_sympy = PermutationGroup(*sympy_perms)
            all_subgroups = _enumerate_all_sympy_subgroups(G_sympy)

            g_order = G_sympy.order()
            seen_fp = set()
            N = 0
            for H in all_subgroups:
                if H.order() == g_order:  # 跳过全群（已单独压入）
                    continue
                if H.order() == 1:  # 跳过恒等群（将单独压入作为格底）
                    continue
                generators = [_sympy_perm_to_tuple(p, self.n) for p in H.generators]
                orbits = self._compute_orbits_by_gen(generators)
                fp = frozenset(frozenset(orbit) for orbit in orbits)
                if fp in seen_fp:
                    continue
                seen_fp.add(fp)
                heapq.heappush(task_heap, (len(orbits), _task_idx, {
                    "name": f"Subgroup |H|={H.order()} (Lattice)",
                    "orbits": orbits,
                    "score": len(orbits)
                }))
                _task_idx += 1
                N += 1

            if self.verbose:
                print(f"  [Lattice] Found {N} distinct subgroup orbit-partitions (|G|={len(self.perms)})")

        except Exception as e:
            print(f"  [Warning] SymPy subgroup enumeration failed: {e}. Falling back to cyclic subgroups.")
            seen_fp = set()
            for p in self.perms:
                orbits = self._compute_orbits_by_gen([p])
                fp = tuple(sorted([len(o) for o in orbits]))
                if fp not in seen_fp:
                    seen_fp.add(fp)
                    heapq.heappush(task_heap, (len(orbits), _task_idx, {
                        "name": f"Cyclic Type {fp[:2]}...",
                        "orbits": orbits,
                        "score": len(orbits)
                    }))
                    _task_idx += 1

        # 3. 恒等群 (Identity) — 子群格底端
        identity_orbits = [[i] for i in range(self.n)]
        heapq.heappush(task_heap, (self.n, _task_idx, {"name": "Identity ({e})", "orbits": identity_orbits, "score": self.n}))
        _task_idx += 1

        if self.verbose:
            print(f"  [Lattice] Generated {len(task_heap)} tasks. Executing full subgroup lattice descent...")

        # 4. 执行：按 orbit count 升序弹出（高对称优先），所有任务运行至穷尽
        initial_active = list(range(1, len(self.perms)))

        while task_heap:
            _, _, task = heapq.heappop(task_heap)
            name = task['name']
            orbits = task['orbits']
            n_orbits = len(orbits)

            if self.verbose:
                print(f"    [Run] {name} (Orbits: {n_orbits})...")

            current_coloring = [None] * self.n
            self._dfs_lattice_pruning(0, orbits, atom_types, current_coloring, initial_active)

    # ---------------- 主入口 ----------------
    def enumerate(self, atom_types: List[str], composition: Dict[str, int] = None):
        self._trie_root = _TrieNode()
        self.classes = {}

        # 统一执行：子群格下降（无截断，保证 Burnside 完备性）
        self.enumerate_subgroups(atom_types)

        # 排序输出
        res = sorted(list(self.classes.values()), key=lambda x: x['stabilizer_size'], reverse=True)
        for i, r in enumerate(res): r['id'] = i

        self._verify_burnside()
        return res

    def _verify_burnside(self):
        lhs = sum(1.0 / d['stabilizer_size'] for d in self.classes.values())
        group_size = len(self.perms)
        if not hasattr(self, '_last_atom_types') or group_size == 0:
            print(f"  [Verification] Sum(1/|Stab|) = {lhs:.4f}")
            return

        atom_types = self._last_atom_types
        # 与 BridgeVerifier.compute_rhs_theoretical 保持完全一致：
        # 有度约束时用 Π |ligands[d]|^count_d，否则用 |atom_types|^n
        if self.topology_requirements and self.user_selected_ligands:
            total_schemes = 1
            for deg, count in self.topology_requirements.items():
                n_choices = len(self.user_selected_ligands.get(deg, atom_types))
                total_schemes *= (n_choices ** count)
        else:
            total_schemes = len(atom_types) ** self.n

        rhs = total_schemes / group_size
        diff = abs(lhs - rhs)
        if diff < 1e-9:
            print(f"  [Verification] [SUCCESS] Sum(1/|Stab|) = {lhs:.6f} = |X|/|G|")
        else:
            coverage = lhs / rhs * 100 if rhs > 0 else 0.0
            print(f"  [Verification] [PARTIAL] Sum(1/|Stab|) = {lhs:.6f}")
            print(f"  [Verification]           |X|/|G|       = {rhs:.6f}")
            print(f"  [Verification]           已覆盖比例     = {coverage:.4f}%")
            print(f"  [Verification]           缺失Burnside量 = {rhs - lhs:.6f}  ≈ 未枚举的低对称轨道数")