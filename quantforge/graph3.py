"""Graph centrality measures: PageRank, degree, closeness, betweenness.

Different notions of how "important" a node is. PageRank ranks by the stationary
distribution of a random surfer (with teleportation); degree centrality counts
neighbors; closeness rewards short distances to everyone; and Brandes' betweenness
counts how often a node sits on shortest paths. These drive link analysis, network
influence, and bottleneck detection. Graphs are ``{node: [neighbors]}`` (or
``{node: {neighbor: weight}}``; weights are ignored except by PageRank's out-degree).
Pure standard library.
"""

from collections import deque


def pagerank(graph, damping=0.85, tol=1e-10, max_iter=1000):
    """PageRank scores by power iteration with teleportation.

    ``graph`` maps each node to its out-neighbors (list or ``{neighbor: weight}``).
    ``damping`` is the follow-a-link probability (``1 - damping`` teleports uniformly).
    Dangling nodes (no out-links) redistribute their mass uniformly. Returns
    ``{node: score}`` summing to 1; converges when the L1 change drops below ``tol``.
    """
    nodes = list(graph)
    n = len(nodes)
    if n == 0:
        raise ValueError("graph must be non-empty")
    if not (0.0 < damping < 1.0):
        raise ValueError("damping must be in (0, 1)")
    out = {u: list(graph[u]) for u in nodes}
    rank = {u: 1.0 / n for u in nodes}
    for _ in range(max_iter):
        dangling = sum(rank[u] for u in nodes if not out[u])
        new = {}
        base = (1.0 - damping) / n + damping * dangling / n
        for v in nodes:
            new[v] = base
        for u in nodes:
            deg = len(out[u])
            if deg:
                share = damping * rank[u] / deg
                for v in out[u]:
                    if v in new:
                        new[v] += share
        diff = sum(abs(new[u] - rank[u]) for u in nodes)
        rank = new
        if diff < tol:
            break
    return rank


def degree_centrality(graph):
    """Degree centrality: each node's neighbor count normalized by ``n - 1``.

    For an undirected graph pass a symmetric adjacency. Returns ``{node: centrality}`` in
    ``[0, 1]`` (a node linked to all others scores 1).
    """
    nodes = list(graph)
    n = len(nodes)
    if n <= 1:
        return {u: 0.0 for u in nodes}
    return {u: len(graph[u]) / (n - 1) for u in nodes}


def closeness_centrality(graph):
    """Closeness centrality: inverse of the mean shortest-path distance (unweighted).

    For each node, ``(reachable) / sum(distances)`` scaled by the fraction reachable
    (Wasserman-Faust), so isolated or unreachable nodes score low. Uses BFS from every
    node. Returns ``{node: centrality}``.
    """
    nodes = list(graph)
    n = len(nodes)
    result = {}
    for s in nodes:
        dist = {s: 0}
        q = deque([s])
        while q:
            u = q.popleft()
            for v in graph.get(u, ()):
                if v not in dist:
                    dist[v] = dist[u] + 1
                    q.append(v)
        total = sum(dist.values())
        reachable = len(dist) - 1
        if total > 0 and n > 1:
            result[s] = (reachable / total) * (reachable / (n - 1))
        else:
            result[s] = 0.0
    return result


def betweenness_centrality(graph, normalized=True):
    """Betweenness centrality by Brandes' algorithm (unweighted, undirected/directed).

    The fraction of shortest paths (between all pairs) that pass through each node.
    ``normalized`` divides by the number of pairs so scores are comparable across graph
    sizes. Returns ``{node: centrality}``.
    """
    nodes = list(graph)
    betw = {u: 0.0 for u in nodes}
    for s in nodes:
        stack = []
        pred = {u: [] for u in nodes}
        sigma = {u: 0.0 for u in nodes}
        sigma[s] = 1.0
        dist = {u: -1 for u in nodes}
        dist[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            stack.append(u)
            for w in graph.get(u, ()):
                if dist[w] < 0:
                    dist[w] = dist[u] + 1
                    q.append(w)
                if dist[w] == dist[u] + 1:
                    sigma[w] += sigma[u]
                    pred[w].append(u)
        delta = {u: 0.0 for u in nodes}
        while stack:
            w = stack.pop()
            for v in pred[w]:
                if sigma[w] > 0:
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            if w != s:
                betw[w] += delta[w]
    n = len(nodes)
    if normalized and n > 2:
        scale = 1.0 / ((n - 1) * (n - 2))
        # undirected counts each pair twice; leave directed-style scaling here and
        # halve below only if the graph is symmetric is caller's concern.
        for u in betw:
            betw[u] *= scale
    return betw
