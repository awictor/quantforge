"""Graph algorithms IV: Bellman-Ford, Floyd-Warshall, A* search.

Shortest paths beyond Dijkstra: Bellman-Ford handles negative edge weights and detects
negative cycles; Floyd-Warshall gives all-pairs shortest distances in ``O(V^3)``; and A*
speeds up point-to-point search with an admissible heuristic. Graphs are
``{node: {neighbor: weight}}``. Pure standard library.
"""

import heapq


def bellman_ford(graph, source):
    """Single-source shortest paths allowing negative edge weights (Bellman-Ford).

    ``graph`` is ``{node: {neighbor: weight}}``; weights may be negative. Returns
    ``(distances, predecessors)`` after ``V - 1`` relaxation rounds. Raises
    ``ValueError`` if a negative-weight cycle is reachable (no finite shortest path).
    """
    if source not in graph:
        raise ValueError("source not in graph")
    nodes = set(graph)
    for u in graph:
        nodes.update(graph[u])
    dist = {n: float("inf") for n in nodes}
    prev = {n: None for n in nodes}
    dist[source] = 0.0
    for _ in range(len(nodes) - 1):
        changed = False
        for u in graph:
            if dist[u] == float("inf"):
                continue
            for v, w in graph[u].items():
                if dist[u] + w < dist[v]:
                    dist[v] = dist[u] + w
                    prev[v] = u
                    changed = True
        if not changed:
            break
    # One more pass: any relaxation means a negative cycle.
    for u in graph:
        if dist[u] == float("inf"):
            continue
        for v, w in graph[u].items():
            if dist[u] + w < dist[v] - 1e-12:
                raise ValueError("graph contains a negative-weight cycle")
    return dist, prev


def floyd_warshall(graph):
    """All-pairs shortest distances by Floyd-Warshall (``O(V^3)``).

    ``graph`` is ``{node: {neighbor: weight}}`` (negative edges allowed, no negative
    cycles). Returns a nested dict ``dist[u][v]`` of shortest path weights (``inf`` if
    unreachable, ``0`` on the diagonal). Raises if a negative cycle is present (a
    diagonal entry goes negative).
    """
    nodes = set(graph)
    for u in graph:
        nodes.update(graph[u])
    nodes = list(nodes)
    dist = {u: {v: float("inf") for v in nodes} for u in nodes}
    for u in nodes:
        dist[u][u] = 0.0
    for u in graph:
        for v, w in graph[u].items():
            if w < dist[u][v]:
                dist[u][v] = w
    for k in nodes:
        dk = dist[k]
        for i in nodes:
            dik = dist[i][k]
            if dik == float("inf"):
                continue
            di = dist[i]
            for j in nodes:
                nd = dik + dk[j]
                if nd < di[j]:
                    di[j] = nd
    for u in nodes:
        if dist[u][u] < -1e-12:
            raise ValueError("graph contains a negative-weight cycle")
    return dist


def a_star(graph, source, target, heuristic):
    """Point-to-point shortest path by A* search with an admissible ``heuristic``.

    ``graph`` is ``{node: {neighbor: weight}}`` with non-negative weights. ``heuristic``
    is a callable ``h(node)`` estimating the remaining cost to ``target``; it must never
    overestimate (admissible) for the result to be optimal. Returns ``(path, cost)``, or
    ``(None, inf)`` if ``target`` is unreachable.
    """
    if source not in graph:
        raise ValueError("source not in graph")
    g = {source: 0.0}
    prev = {source: None}
    pq = [(heuristic(source), source)]
    closed = set()
    while pq:
        _, u = heapq.heappop(pq)
        if u == target:
            path = []
            node = u
            while node is not None:
                path.append(node)
                node = prev[node]
            path.reverse()
            return path, g[u]
        if u in closed:
            continue
        closed.add(u)
        for v, w in graph.get(u, {}).items():
            if w < 0:
                raise ValueError("A* requires non-negative weights")
            ng = g[u] + w
            if v not in g or ng < g[v]:
                g[v] = ng
                prev[v] = u
                heapq.heappush(pq, (ng + heuristic(v), v))
    return None, float("inf")
