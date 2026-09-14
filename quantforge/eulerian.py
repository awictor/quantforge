"""Eulerian paths and circuits (Hierholzer's algorithm).

An Eulerian *circuit* traverses every edge of a graph exactly once and returns to its
start; an Eulerian *path* does the same without the return requirement. They exist only
under degree conditions -- for an undirected connected graph, a circuit needs every degree
even, a path needs exactly zero or two odd-degree vertices; for a directed graph, in-degree
must equal out-degree everywhere (circuit) or differ by one at a single source/sink pair
(path). Hierholzer's algorithm builds the trail in ``O(E)`` by splicing cycles. Pure
standard library.
"""

from collections import defaultdict


def _degrees_undirected(graph):
    deg = defaultdict(int)
    edge_count = 0
    for u, nbrs in graph.items():
        deg[u] += len(nbrs)
        for v in nbrs:
            deg[v]  # ensure v is present
            edge_count += 1
    return deg, edge_count // 2


def _connected_ignoring_isolated(vertices, adj):
    # check all vertices with degree > 0 are in one connected component
    active = [v for v in vertices if adj[v]]
    if not active:
        return True
    seen = {active[0]}
    stack = [active[0]]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return all(v in seen for v in active)


def has_eulerian_path(graph, directed=False):
    """True if ``graph`` has an Eulerian path (trail using every edge once)."""
    return _euler_kind(graph, directed) in ("path", "circuit")


def has_eulerian_circuit(graph, directed=False):
    """True if ``graph`` has an Eulerian circuit (closed Eulerian path)."""
    return _euler_kind(graph, directed) == "circuit"


def _all_vertices(graph):
    verts = set(graph)
    for u in graph:
        verts.update(graph[u])
    return verts


def _euler_kind(graph, directed):
    # returns "circuit", "path", or None
    verts = _all_vertices(graph)
    adj = {v: list(graph.get(v, [])) for v in verts}
    total_edges = sum(len(a) for a in adj.values())
    if total_edges == 0:
        return "circuit"          # empty graph: trivially a (degenerate) circuit

    if directed:
        indeg = defaultdict(int)
        outdeg = defaultdict(int)
        for u in adj:
            outdeg[u] += len(adj[u])
            for v in adj[u]:
                indeg[v] += 1
        start_extra = end_extra = 0
        for v in verts:
            diff = outdeg[v] - indeg[v]
            if diff == 1:
                start_extra += 1
            elif diff == -1:
                end_extra += 1
            elif diff != 0:
                return None
        # weak connectivity of edge-bearing vertices
        undirected = {v: [] for v in verts}
        for u in adj:
            for v in adj[u]:
                undirected[u].append(v)
                undirected[v].append(u)
        if not _connected_ignoring_isolated(verts, undirected):
            return None
        if start_extra == 0 and end_extra == 0:
            return "circuit"
        if start_extra == 1 and end_extra == 1:
            return "path"
        return None
    else:
        deg = {v: len(adj[v]) for v in verts}
        # each undirected edge appears once here (graph stores one direction); but to count
        # true degree we need both endpoints. Rebuild symmetric degree:
        sym = defaultdict(int)
        for u in adj:
            for v in adj[u]:
                sym[u] += 1
                sym[v] += 1
        odd = [v for v in verts if sym[v] % 2 == 1]
        undirected = {v: [] for v in verts}
        for u in adj:
            for v in adj[u]:
                undirected[u].append(v)
                undirected[v].append(u)
        if not _connected_ignoring_isolated(verts, undirected):
            return None
        if len(odd) == 0:
            return "circuit"
        if len(odd) == 2:
            return "path"
        return None


def eulerian_path(graph, directed=False, start=None):
    """Return a list of vertices tracing an Eulerian path, or ``None`` if none exists.

    The returned list has ``E + 1`` entries (each consecutive pair is a traversed edge).
    For a circuit the first and last vertices coincide. ``graph`` is ``{u: [v, ...]}``; for
    an undirected graph list each edge once (both directions are inferred).
    """
    kind = _euler_kind(graph, directed)
    if kind is None:
        return None

    verts = _all_vertices(graph)
    if directed:
        adj = {v: list(reversed(graph.get(v, []))) for v in verts}  # pop() takes last
    else:
        adj = {v: [] for v in verts}
        # build a multigraph with per-edge "used" flags shared between endpoints
        edge_id = 0
        used = {}
        for u in graph:
            for v in graph[u]:
                adj[u].append((v, edge_id))
                adj[v].append((u, edge_id))
                used[edge_id] = False
                edge_id += 1

    total_edges = sum(len(graph.get(v, [])) for v in verts)
    if total_edges == 0:
        # degenerate: return the lone start vertex if any
        s = start if start is not None else (next(iter(verts)) if verts else None)
        return [s] if s is not None else []

    # choose start
    if start is None:
        start = _default_start(graph, directed, kind)

    stack = [start]
    trail = []
    if directed:
        while stack:
            u = stack[-1]
            if adj[u]:
                stack.append(adj[u].pop())
            else:
                trail.append(stack.pop())
    else:
        ptr = {v: 0 for v in verts}
        while stack:
            u = stack[-1]
            advanced = False
            lst = adj[u]
            while ptr[u] < len(lst):
                v, eid = lst[ptr[u]]
                ptr[u] += 1
                if not used[eid]:
                    used[eid] = True
                    stack.append(v)
                    advanced = True
                    break
            if not advanced:
                trail.append(stack.pop())
    trail.reverse()
    return trail


def _default_start(graph, directed, kind):
    verts = _all_vertices(graph)
    if directed:
        indeg = defaultdict(int)
        outdeg = defaultdict(int)
        for u in graph:
            for v in graph[u]:
                outdeg[u] += 1
                indeg[v] += 1
        if kind == "path":
            for v in verts:
                if outdeg[v] - indeg[v] == 1:
                    return v
        for v in verts:
            if outdeg[v] > 0:
                return v
    else:
        sym = defaultdict(int)
        for u in graph:
            for v in graph[u]:
                sym[u] += 1
                sym[v] += 1
        if kind == "path":
            for v in verts:
                if sym[v] % 2 == 1:
                    return v
        for v in verts:
            if sym[v] > 0:
                return v
    return next(iter(verts)) if verts else None
