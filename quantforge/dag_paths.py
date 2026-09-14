"""Longest/shortest paths and reachability on a directed acyclic graph.

On a DAG the vertices can be linearized (topological order), which turns path optimization
into a single ``O(V + E)`` sweep -- and unlike Dijkstra it handles *negative* edge weights,
and unlike Bellman-Ford it needs no relaxation rounds. The *longest* path (critical path) is
the same sweep maximizing instead of minimizing, the backbone of task scheduling (CPM).
``transitive_closure`` gives all-pairs reachability. Raises if the graph has a cycle. Pure
standard library.
"""


def _topological_order(graph):
    indeg = {u: 0 for u in graph}
    for u in graph:
        for v, _w in graph[u]:
            if v not in indeg:
                raise ValueError("edge to unknown node %r" % (v,))
            indeg[v] += 1
    queue = [u for u in graph if indeg[u] == 0]
    order = []
    i = 0
    while i < len(queue):
        u = queue[i]
        i += 1
        order.append(u)
        for v, _w in graph[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                queue.append(v)
    if len(order) != len(graph):
        raise ValueError("graph has a cycle; not a DAG")
    return order


def _extremal_path(graph, source, longest):
    if source not in graph:
        raise ValueError("source not in graph")
    order = _topological_order(graph)
    INF = float("-inf") if longest else float("inf")
    dist = {u: INF for u in graph}
    dist[source] = 0.0
    prev = {u: None for u in graph}
    started = False
    for u in order:
        if u == source:
            started = True
        if not started or dist[u] == INF:
            continue
        for v, w in graph[u]:
            nd = dist[u] + w
            better = nd > dist[v] if longest else nd < dist[v]
            if better:
                dist[v] = nd
                prev[v] = u
    return dist, prev


def dag_shortest_path(graph, source, target=None):
    """Shortest-path distances from ``source`` on a DAG (negative weights allowed).

    ``graph`` is ``{node: [(neighbor, weight), ...]}``. Returns the distance dict; if
    ``target`` is given, returns ``(distance, path)`` for that target instead
    (``distance`` is ``inf`` and ``path`` empty if unreachable).
    """
    dist, prev = _extremal_path(graph, source, longest=False)
    if target is None:
        return dist
    return _reconstruct(dist, prev, source, target, float("inf"))


def dag_longest_path(graph, source, target=None):
    """Longest-path distances from ``source`` on a DAG (the critical path).

    Same arguments as :func:`dag_shortest_path`. Unreachable nodes have distance
    ``-inf``.
    """
    dist, prev = _extremal_path(graph, source, longest=True)
    if target is None:
        return dist
    return _reconstruct(dist, prev, source, target, float("-inf"))


def _reconstruct(dist, prev, source, target, unreachable):
    if target not in dist or dist[target] == unreachable:
        return unreachable, []
    path = []
    u = target
    while u is not None:
        path.append(u)
        if u == source:
            break
        u = prev[u]
    path.reverse()
    return dist[target], path


def transitive_closure(graph):
    """Return ``{node: set(reachable nodes)}`` (excluding the node itself unless it loops).

    ``graph`` may be ``{node: [(neighbor, weight), ...]}`` or ``{node: [neighbor, ...]}``.
    Computed by a DFS from each node; works on any directed graph (cycles allowed).
    """
    adj = {}
    for u in graph:
        succ = []
        for e in graph[u]:
            succ.append(e[0] if isinstance(e, (tuple, list)) else e)
        adj[u] = succ

    def reach(start):
        seen = set()
        stack = list(adj.get(start, []))
        while stack:
            v = stack.pop()
            if v not in seen:
                seen.add(v)
                stack.extend(adj.get(v, []))
        return seen

    return {u: reach(u) for u in adj}
