"""Graph algorithms II: union-find, minimum spanning tree, maximum flow.

Weighted-graph structure problems: a union-find (disjoint-set) forest with path
compression and union by rank, Kruskal's minimum spanning tree built on it, and
Edmonds-Karp maximum flow (BFS-augmenting Ford-Fulkerson). These solve network design,
clustering, and capacity/bottleneck questions. Pure standard library.
"""

from collections import deque


class UnionFind:
    """Disjoint-set forest with path compression and union by rank.

    ``find(x)`` returns the representative of ``x``'s set; ``union(a, b)`` merges two
    sets and returns whether they were previously disjoint; ``connected(a, b)`` tests
    membership. Elements are created on first reference. Near-constant amortized cost.
    """

    def __init__(self):
        self.parent = {}
        self.rank = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            return x
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        # Path compression.
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True

    def connected(self, a, b):
        return self.find(a) == self.find(b)


def minimum_spanning_tree(nodes, edges):
    """Minimum spanning tree (or forest) by Kruskal's algorithm.

    ``nodes`` is an iterable of vertices; ``edges`` a list of ``(u, v, weight)``. Returns
    ``(tree_edges, total_weight)`` where ``tree_edges`` is the chosen subset (sorted by
    weight then endpoints). For a disconnected graph this is the minimum spanning forest.
    """
    uf = UnionFind()
    for n in nodes:
        uf.find(n)
    tree = []
    total = 0.0
    for u, v, w in sorted(edges, key=lambda e: (e[2], e[0], e[1])):
        if uf.union(u, v):
            tree.append((u, v, w))
            total += w
    return tree, total


def max_flow(graph, source, sink):
    """Maximum flow from ``source`` to ``sink`` (Edmonds-Karp).

    ``graph`` is ``{node: {neighbor: capacity}}`` with capacities ``>= 0``. Returns the
    maximum total flow value. Uses BFS to find shortest augmenting paths in the residual
    network, so it runs in ``O(V E^2)`` and terminates on rational capacities.
    """
    if source == sink:
        raise ValueError("source and sink must differ")
    # Residual capacities; ensure reverse arcs exist.
    residual = {}
    for u in graph:
        residual.setdefault(u, {})
        for v, cap in graph[u].items():
            if cap < 0:
                raise ValueError("capacities must be non-negative")
            residual[u][v] = residual[u].get(v, 0) + cap
            residual.setdefault(v, {})
            residual[v].setdefault(u, 0)
    flow = 0.0
    while True:
        # BFS for an augmenting path, tracking parents.
        parent = {source: None}
        q = deque([source])
        while q and sink not in parent:
            u = q.popleft()
            for v, cap in residual[u].items():
                if v not in parent and cap > 0:
                    parent[v] = u
                    q.append(v)
        if sink not in parent:
            break
        # Bottleneck along the path.
        bottleneck = float("inf")
        v = sink
        while parent[v] is not None:
            u = parent[v]
            bottleneck = min(bottleneck, residual[u][v])
            v = u
        # Augment.
        v = sink
        while parent[v] is not None:
            u = parent[v]
            residual[u][v] -= bottleneck
            residual[v][u] += bottleneck
            v = u
        flow += bottleneck
    return flow
