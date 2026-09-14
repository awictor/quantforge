"""Minimum-cost maximum-flow (successive shortest augmenting paths).

Given a network with per-edge capacity *and* cost, this finds a maximum-flow from source to
sink that, among all maximum flows, has the least total cost. It augments along a
shortest-cost path each round (found by Bellman-Ford/SPFA, so negative-cost edges are fine),
pushing as much flow as the path's residual allows, until the sink is unreachable. Returns
the total flow and its cost. The plain maximum flow is the special case of all-equal costs.
Pure standard library.
"""

from collections import deque


class MinCostMaxFlow:
    """Min-cost max-flow on an integer-capacity network built edge by edge.

    ``add_edge(u, v, capacity, cost)`` adds a directed edge (and its residual). ``solve(s,
    t)`` returns ``(max_flow, min_cost)``: the maximum flow value and the least total cost
    achieving it. Nodes are any hashable labels.
    """

    __slots__ = ("_graph", "_edges")

    def __init__(self):
        self._graph = {}          # node -> list of edge indices
        self._edges = []          # [to, capacity, cost, flow]

    def _node(self, u):
        if u not in self._graph:
            self._graph[u] = []

    def add_edge(self, u, v, capacity, cost):
        """Add a directed edge ``u -> v`` with the given non-negative ``capacity`` and ``cost``."""
        if capacity < 0:
            raise ValueError("capacity must be non-negative")
        self._node(u)
        self._node(v)
        self._graph[u].append(len(self._edges))
        self._edges.append([v, capacity, cost, 0])
        self._graph[v].append(len(self._edges))
        self._edges.append([u, 0, -cost, 0])   # residual edge
        return self

    def solve(self, source, sink):
        """Return ``(max_flow, min_cost)`` from ``source`` to ``sink``."""
        if source == sink:
            raise ValueError("source and sink must differ")
        if source not in self._graph or sink not in self._graph:
            return 0, 0            # a node with no incident edges carries no flow
        total_flow = 0
        total_cost = 0
        INF = float("inf")
        while True:
            # SPFA: shortest-cost path in the residual graph
            dist = {u: INF for u in self._graph}
            dist[source] = 0
            in_queue = {u: False for u in self._graph}
            prev_edge = {u: -1 for u in self._graph}
            q = deque([source])
            in_queue[source] = True
            while q:
                u = q.popleft()
                in_queue[u] = False
                du = dist[u]
                for ei in self._graph[u]:
                    to, cap, cost, flow = self._edges[ei]
                    if cap - flow > 0 and du + cost < dist[to]:
                        dist[to] = du + cost
                        prev_edge[to] = ei
                        if not in_queue[to]:
                            q.append(to)
                            in_queue[to] = True
            if dist[sink] == INF:
                break                      # sink unreachable: done
            # bottleneck along the found path
            push = INF
            v = sink
            while v != source:
                ei = prev_edge[v]
                to, cap, cost, flow = self._edges[ei]
                push = min(push, cap - flow)
                v = self._edges[ei ^ 1][0]
            # apply the augmentation
            v = sink
            while v != source:
                ei = prev_edge[v]
                self._edges[ei][3] += push
                self._edges[ei ^ 1][3] -= push
                v = self._edges[ei ^ 1][0]
            total_flow += push
            total_cost += push * dist[sink]
        return total_flow, total_cost


def min_cost_max_flow(edges, source, sink):
    """Convenience: build the network from ``edges`` = ``[(u, v, capacity, cost), ...]``.

    Returns ``(max_flow, min_cost)``.
    """
    mcmf = MinCostMaxFlow()
    for u, v, cap, cost in edges:
        mcmf.add_edge(u, v, cap, cost)
    return mcmf.solve(source, sink)
