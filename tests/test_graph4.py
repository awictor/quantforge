"""Graph algorithms IV: Bellman-Ford, Floyd-Warshall, A* search."""

import random

import pytest

from quantforge import bellman_ford, floyd_warshall, a_star
from quantforge import dijkstra


def test_bellman_ford_negative_edge():
    g = {'A': {'B': 4, 'C': 5}, 'B': {'C': -3}, 'C': {'D': 2}, 'D': {}}
    dist, _ = bellman_ford(g, 'A')
    assert dist == {'A': 0, 'B': 4, 'C': 1, 'D': 3}


def test_bellman_ford_detects_negative_cycle():
    with pytest.raises(ValueError):
        bellman_ford({'A': {'B': 1}, 'B': {'C': -3}, 'C': {'A': 1}}, 'A')


def test_bellman_ford_matches_dijkstra():
    rng = random.Random(1)
    for _ in range(100):
        nodes = list(range(6))
        graph = {n: {} for n in nodes}
        for u in nodes:
            for v in nodes:
                if u != v and rng.random() < 0.4:
                    graph[u][v] = rng.randint(1, 9)
        bf, _ = bellman_ford(graph, 0)
        dj, _ = dijkstra(graph, 0)
        for n in nodes:
            assert bf[n] == dj[n]


def test_floyd_warshall_matches_dijkstra():
    rng = random.Random(2)
    for _ in range(50):
        nodes = list(range(6))
        graph = {n: {} for n in nodes}
        for u in nodes:
            for v in nodes:
                if u != v and rng.random() < 0.4:
                    graph[u][v] = rng.randint(1, 9)
        fw = floyd_warshall(graph)
        for s in nodes:
            dj, _ = dijkstra(graph, s)
            for t in nodes:
                a, b = fw[s][t], dj[t]
                assert a == b or abs(a - b) < 1e-9


def _grid(n):
    g = {}
    for i in range(n):
        for j in range(n):
            nb = {}
            for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                ni, nj = i + di, j + dj
                if 0 <= ni < n and 0 <= nj < n:
                    nb[(ni, nj)] = 1
            g[(i, j)] = nb
    return g


def test_astar_optimal_on_grid():
    G = _grid(5)
    h = lambda p: abs(p[0] - 4) + abs(p[1] - 4)
    path, cost = a_star(G, (0, 0), (4, 4), h)
    dj, _ = dijkstra(G, (0, 0))
    assert cost == dj[(4, 4)] == 8
    assert path[0] == (0, 0) and path[-1] == (4, 4)


def test_astar_zero_heuristic_and_unreachable():
    g = {'A': {'B': 1, 'C': 4}, 'B': {'C': 2, 'D': 5}, 'C': {'D': 1}, 'D': {}}
    assert a_star(g, 'A', 'D', lambda n: 0) == (['A', 'B', 'C', 'D'], 4.0)
    assert a_star({'A': {'B': 1}, 'B': {}, 'C': {}}, 'A', 'C', lambda n: 0) == \
        (None, float('inf'))


def test_validation():
    with pytest.raises(ValueError):
        bellman_ford({'A': {}}, 'Z')
    with pytest.raises(ValueError):
        a_star({'A': {'B': -1}, 'B': {}}, 'A', 'B', lambda n: 0)
