"""Coverage certificates for completed DFS subtrees of interrupted tasks.

Nodes are constructed only from an interrupted DFS frontier. A terminal means
that the canonical search under that task's partition and prefix is exhausted.
It says nothing about an arbitrary matching prefix in a different task: the
query also proves containment in the source fixed-point domain and label pools.
"""


class _PrefixNode:
    __slots__ = ("children", "closed")

    def __init__(self):
        self.children = {}
        self.closed = False


def _vertex_domains(orbits, domains):
    n = sum(map(len, orbits))
    by_node = [None] * n
    for orbit, domain in zip(orbits, domains):
        choices = frozenset(domain)
        for node in orbit:
            by_node[node] = choices
    return tuple(by_node)


class CompletedPrefixCache:
    def __init__(self, max_nodes=50000):
        if max_nodes < 0:
            raise ValueError("prefix_trie_max_nodes must be non-negative")
        self.max_nodes = max_nodes
        self.forests = []
        self.stats = dict(queries=0, hits=0, early_hits=0, closed_prefixes=0,
                          nodes_created=0, live_nodes=0, peak_nodes=0,
                          budget_skips=0, snapshots=0, forests=0)

    def _new(self):
        if self.stats["live_nodes"] >= self.max_nodes:
            self.stats["budget_skips"] += 1
            return None
        self.stats["nodes_created"] += 1
        self.stats["live_nodes"] += 1
        self.stats["peak_nodes"] = max(self.stats["peak_nodes"], self.stats["live_nodes"])
        return _PrefixNode()

    def _child(self, node, value):
        child = node.children.get(value)
        if child is None:
            child = self._new()
            if child is not None:
                node.children[value] = child
        return child

    @staticmethod
    def _count(node):
        stack = [node]
        count = 0
        while stack:
            current = stack.pop()
            count += 1
            stack.extend(current.children.values())
        return count

    def capture(self, orbits, domains, values, progress, depth):
        """Save completed siblings of the interrupted path, never the path itself.

        progress[k] is the number of exhausted labels at the active depth k;
        values[:depth] is the interrupted path. Deeper/stale frames are ignored.
        """
        self.stats["snapshots"] += 1
        root = self._new()
        if root is None:
            return
        node = root
        closed = 0
        for level in range(min(depth + 1, len(orbits))):
            for value in domains[level][:progress[level]]:
                child = self._child(node, value)
                if child is not None and not child.closed:
                    self.stats["live_nodes"] -= sum(self._count(x) for x in child.children.values())
                    child.children.clear()
                    child.closed = True
                    closed += 1
            if level >= depth or node.closed:
                break
            child = self._child(node, values[level])
            if child is None or child.closed:
                break
            node = child
        if closed:
            self.forests.append((tuple(tuple(orbit) for orbit in orbits),
                                 _vertex_domains(orbits, domains), root))
            self.stats["closed_prefixes"] += closed
            self.stats["forests"] = len(self.forests)
        else:
            self.stats["live_nodes"] -= self._count(root)

    def task_query(self, orbits, domains):
        current_domains = _vertex_domains(orbits, domains)
        eligible = [(partition, root) for partition, source_domains, root in self.forests
                    if len(current_domains) == len(source_domains)
                    and all(current <= source for current, source in zip(current_domains, source_domains))]
        return _TaskQuery(self, orbits, domains, eligible)


class _TaskQuery:
    def __init__(self, cache, orbits, domains, forests):
        self.cache = cache
        self.forests = forests
        self.depth_limit = len(orbits)
        self.at_depth = [[] for _ in range(len(orbits) + 1)]
        self.expressions = [None] * sum(map(len, orbits))
        for index, (orbit, domain) in enumerate(zip(orbits, domains)):
            # Nonnegative expressions are fixed values; negative expressions
            # identify an unassigned current-task variable shared by its orbit.
            expression = domain[0] if len(domain) == 1 else -index-1
            for node in orbit:
                self.expressions[node] = expression
        for partition, root in forests:
            equality_depth = 0
            impossible = False
            for orbit in partition:
                terms = {self.expressions[node] for node in orbit}
                if len(terms) <= 1:
                    continue
                if sum(term >= 0 for term in terms) > 1:
                    impossible = True
                    break
                equality_depth = max(equality_depth, max((-term for term in terms if term < 0), default=0))
            if impossible:
                continue
            # A closed prefix can first become decidable at one depth. A
            # failed equality or label match stays failed in its descendants.
            trigger_depths = set()
            stack = [(root, 0, equality_depth)]
            while stack:
                node, level, required = stack.pop()
                if node.closed:
                    trigger_depths.add(required)
                elif node.children:
                    term = self.expressions[partition[level][0]]
                    required = max(required, -term if term < 0 else 0)
                    stack.extend((child, level + 1, required) for child in node.children.values())
            for depth in trigger_depths:
                self.at_depth[depth].append((partition, root))

    def covers(self, coloring, depth):
        if not self.at_depth[depth]:
            return False
        self.cache.stats["queries"] += 1
        expressions = [value if value is not None else unknown
                       for value, unknown in zip(coloring, self.expressions)]
        for partition, root in self.at_depth[depth]:
            # Every completion must satisfy all source-task equalities.
            if any(any(expressions[node] != expressions[orbit[0]] for node in orbit[1:])
                   for orbit in partition):
                continue
            current = root
            for orbit in partition:
                value = expressions[orbit[0]]
                if current.closed or value < 0:
                    break
                current = current.children.get(value)
                if current is None:
                    break
            if current is not None and current.closed:
                self.cache.stats["hits"] += 1
                self.cache.stats["early_hits"] += depth < self.depth_limit
                return True
        return False
