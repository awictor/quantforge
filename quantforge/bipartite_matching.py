"""Maximum bipartite matching (Hopcroft-Karp) and related structure.

A matching pairs left vertices with right vertices so no vertex is used twice; a *maximum*
matching pairs as many as possible. Hopcroft-Karp finds one in ``O(E * sqrt(V))`` by
repeatedly finding a maximal set of shortest augmenting paths (a BFS to layer the graph,
then DFS to augment along disjoint paths). By Konig's theorem the maximum matching size
equals the minimum vertex cover in a bipartite graph, and the actual cover is recovered
from the alternating-path reachability of unmatched left vertices. Pure standard library.
"""

from collections import deque


def maximum_bipartite_matching(adjacency):
    """Return a maximum matching as a ``{left: right}`` dict (Hopcroft-Karp).

    ``adjacency`` maps each left vertex to an iterable of the right vertices it can pair
    with. Left and right vertex labels live in separate namespaces (they may overlap in
    value without conflict). Only left vertices present as keys are matched.
    """
    INF = float("inf")
    left_nodes = list(adjacency)
    adj = {u: list(vs) for u, vs in adjacency.items()}
    match_left = {}          # left -> right
    match_right = {}         # right -> left
    dist = {}

    def bfs():
        queue = deque()
        for u in left_nodes:
            if u not in match_left:
                dist[u] = 0
                queue.append(u)
            else:
                dist[u] = INF
        found = False
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                w = match_right.get(v)
                if w is None:
                    found = True            # reached a free right vertex: augmenting path exists
                elif dist[w] == INF:
                    dist[w] = dist[u] + 1
                    queue.append(w)
        return found

    def dfs(u):
        for v in adj[u]:
            w = match_right.get(v)
            if w is None or (dist[w] == dist[u] + 1 and dfs(w)):
                match_left[u] = v
                match_right[v] = u
                return True
        dist[u] = INF
        return False

    while bfs():
        for u in left_nodes:
            if u not in match_left:
                dfs(u)
    return dict(match_left)


def maximum_matching_size(adjacency):
    """Size of a maximum bipartite matching."""
    return len(maximum_bipartite_matching(adjacency))


def minimum_vertex_cover(adjacency):
    """Return a minimum vertex cover ``(left_set, right_set)`` via Konig's theorem.

    Its total size equals the maximum matching size. Computed from the maximum matching by
    marking left vertices reachable by alternating paths from unmatched left vertices: the
    cover is the unmarked left vertices plus the marked right vertices.
    """
    adj = {u: list(vs) for u, vs in adjacency.items()}
    match_left = maximum_bipartite_matching(adjacency)
    match_right = {v: u for u, v in match_left.items()}

    # alternating-path DFS from every unmatched left vertex
    visited_left = set()
    visited_right = set()

    def visit(u):
        visited_left.add(u)
        for v in adj[u]:
            if v not in visited_right:
                visited_right.add(v)
                w = match_right.get(v)
                if w is not None and w not in visited_left:
                    visit(w)

    for u in adj:
        if u not in match_left:
            visit(u)

    left_cover = {u for u in adj if u not in visited_left}
    right_cover = set(visited_right)
    return left_cover, right_cover
