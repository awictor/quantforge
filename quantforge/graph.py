"""Graph algorithms: shortest paths, traversal, components, topological order.

The core graph routines on a dict adjacency list -- Dijkstra's shortest paths (with
non-negative weights), breadth-first traversal and unweighted shortest hops, connected
components of an undirected graph, and Kahn's topological sort with cycle detection.
Graphs are ``{node: {neighbor: weight}}`` (or ``{node: [neighbors]}`` for unweighted
helpers). Pure standard library.
"""

import heapq
from collections import deque


def dijkstra(graph, source):
    """Shortest-path distances from ``source`` over non-negative edge weights.

    ``graph`` is ``{node: {neighbor: weight}}`` with weights ``>= 0``. Returns
    ``(distances, predecessors)``: ``distances[node]`` is the shortest total weight from
    ``source`` (``inf`` if unreachable) and ``predecessors[node]`` the previous node on a
    shortest path (``None`` at the source or if unreachable). Raises on a negative weight.
    """
    if source not in graph:
        raise ValueError("source not in graph")
    dist = {node: float("inf") for node in graph}
    prev = {node: None for node in graph}
    dist[source] = 0.0
    pq = [(0.0, source)]
    visited = set()
    while pq:
        d, u = heapq.heappop(pq)
        if u in visited:
            continue
        visited.add(u)
        for v, w in graph.get(u, {}).items():
            if w < 0:
                raise ValueError("dijkstra requires non-negative weights")
            if v not in dist:
                dist[v] = float("inf")
                prev[v] = None
            nd = d + w
            if nd < dist[v]:
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
    return dist, prev


def shortest_path(graph, source, target):
    """Shortest path ``[source, ..., target]`` and its total weight, via Dijkstra.

    Returns ``(path, distance)``; ``(None, inf)`` if ``target`` is unreachable.
    """
    dist, prev = dijkstra(graph, source)
    if target not in dist or dist[target] == float("inf"):
        return None, float("inf")
    path = []
    node = target
    while node is not None:
        path.append(node)
        node = prev[node]
    path.reverse()
    return path, dist[target]


def bfs(graph, source):
    """Breadth-first hop distances from ``source`` (unweighted).

    ``graph`` maps each node to an iterable of neighbors. Returns ``{node: hops}`` for
    every node reachable from ``source`` (``0`` at the source), in increasing hop order.
    """
    if source not in graph:
        raise ValueError("source not in graph")
    hops = {source: 0}
    q = deque([source])
    while q:
        u = q.popleft()
        for v in graph.get(u, ()):
            if v not in hops:
                hops[v] = hops[u] + 1
                q.append(v)
    return hops


def connected_components(graph):
    """Connected components of an undirected graph as a list of node sets.

    ``graph`` maps each node to its neighbors; edges are treated as undirected (a node
    listed as a neighbor is joined even if the reverse edge is absent). Returns the
    components as a list of sorted lists.
    """
    # Build a symmetric adjacency so a one-directional listing still connects.
    adj = {node: set(graph.get(node, ())) for node in graph}
    for u in list(adj):
        for v in adj[u]:
            adj.setdefault(v, set()).add(u)
    seen = set()
    components = []
    for start in adj:
        if start in seen:
            continue
        stack = [start]
        comp = []
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            comp.append(u)
            stack.extend(adj[u] - seen)
        components.append(sorted(comp))
    return components


def topological_sort(graph):
    """Topological ordering of a DAG (Kahn's algorithm); raises if the graph has a cycle.

    ``graph`` is ``{node: [successors]}`` (or ``{node: {succ: weight}}``). Returns a list
    of nodes such that every edge points forward. Raises ``ValueError`` if a cycle makes
    a valid ordering impossible.
    """
    nodes = set(graph)
    for u in graph:
        nodes.update(graph[u])
    indeg = {n: 0 for n in nodes}
    for u in graph:
        for v in graph[u]:
            indeg[v] += 1
    queue = deque(sorted(n for n in nodes if indeg[n] == 0))
    order = []
    while queue:
        u = queue.popleft()
        order.append(u)
        for v in sorted(graph.get(u, ())):
            indeg[v] -= 1
            if indeg[v] == 0:
                queue.append(v)
    if len(order) != len(nodes):
        raise ValueError("graph has a cycle; no topological order exists")
    return order
