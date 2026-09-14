"""Graph algorithms: shortest paths, traversal, components, topological order."""

import random

import pytest

from quantforge import (
    dijkstra,
    shortest_path,
    bfs,
    connected_components,
    topological_sort,
)


def test_dijkstra_known_distances():
    g = {'A': {'B': 1, 'C': 4}, 'B': {'C': 2, 'D': 5}, 'C': {'D': 1}, 'D': {}}
    dist, _ = dijkstra(g, 'A')
    assert dist == {'A': 0, 'B': 1, 'C': 3, 'D': 4}


def test_shortest_path_reconstruction():
    g = {'A': {'B': 1, 'C': 4}, 'B': {'C': 2, 'D': 5}, 'C': {'D': 1}, 'D': {}}
    path, d = shortest_path(g, 'A', 'D')
    assert path == ['A', 'B', 'C', 'D'] and d == 4
    assert shortest_path({'A': {'B': 1}, 'B': {}, 'C': {}}, 'A', 'C') == (None, float('inf'))


def test_dijkstra_matches_bellman_ford():
    rng = random.Random(1)

    def bellman(graph, src):
        dist = {n: float('inf') for n in graph}
        dist[src] = 0
        for _ in range(len(graph) - 1):
            for u in graph:
                for v, w in graph[u].items():
                    if dist[u] + w < dist[v]:
                        dist[v] = dist[u] + w
        return dist

    for _ in range(100):
        nodes = list(range(6))
        graph = {n: {} for n in nodes}
        for u in nodes:
            for v in nodes:
                if u != v and rng.random() < 0.4:
                    graph[u][v] = rng.randint(1, 9)
        dd, _ = dijkstra(graph, 0)
        assert dd == bellman(graph, 0)


def test_bfs_hops():
    gb = {'A': ['B', 'C'], 'B': ['D'], 'C': ['D'], 'D': ['E'], 'E': []}
    assert bfs(gb, 'A') == {'A': 0, 'B': 1, 'C': 1, 'D': 2, 'E': 3}


def test_connected_components():
    gc = {1: [2], 2: [1], 3: [4], 4: [], 5: []}
    assert connected_components(gc) == [[1, 2], [3, 4], [5]]


def test_topological_sort_is_valid():
    dag = {'shirt': ['tie', 'belt'], 'tie': ['jacket'], 'belt': ['jacket'], 'jacket': []}
    order = topological_sort(dag)
    pos = {n: i for i, n in enumerate(order)}
    for u in dag:
        for v in dag[u]:
            assert pos[u] < pos[v]


def test_topological_sort_detects_cycle():
    with pytest.raises(ValueError):
        topological_sort({'a': ['b'], 'b': ['c'], 'c': ['a']})


def test_validation():
    with pytest.raises(ValueError):
        dijkstra({'A': {}}, 'Z')                       # source not in graph
    with pytest.raises(ValueError):
        dijkstra({'A': {'B': -1}, 'B': {}}, 'A')       # negative weight
    with pytest.raises(ValueError):
        bfs({'A': []}, 'Z')
