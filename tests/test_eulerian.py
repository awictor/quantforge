"""Tests for Eulerian path/circuit (Hierholzer), cross-checked by trail verification."""

import random
from collections import Counter

from quantforge.eulerian import (
    eulerian_path,
    has_eulerian_path,
    has_eulerian_circuit,
)


def _edges_undirected(graph):
    c = Counter()
    for u in graph:
        for v in graph[u]:
            c[frozenset((u, v)) if u != v else (u, "self")] += 1
    return c


def _trail_edges_undirected(trail):
    c = Counter()
    for a, b in zip(trail, trail[1:]):
        c[frozenset((a, b)) if a != b else (a, "self")] += 1
    return c


def _edges_directed(graph):
    c = Counter()
    for u in graph:
        for v in graph[u]:
            c[(u, v)] += 1
    return c


def _trail_edges_directed(trail):
    return Counter((a, b) for a, b in zip(trail, trail[1:]))


def test_fuzz_undirected():
    rng = random.Random(91)
    saw_path = saw_none = 0
    for _ in range(3000):
        n = rng.randint(1, 7)
        graph = {i: [] for i in range(n)}
        for _ in range(rng.randint(0, 10)):
            u = rng.randint(0, n - 1)
            v = rng.randint(0, n - 1)
            if u != v:
                graph[u].append(v)
        trail = eulerian_path(graph, directed=False)
        total = sum(len(graph[u]) for u in graph)
        if has_eulerian_path(graph, False):
            assert trail is not None
            assert len(trail) == total + 1
            assert _edges_undirected(graph) == _trail_edges_undirected(trail)
            if has_eulerian_circuit(graph, False) and total > 0:
                assert trail[0] == trail[-1]
            saw_path += 1
        else:
            assert trail is None
            saw_none += 1
    assert saw_path > 0 and saw_none > 0


def test_fuzz_directed():
    rng = random.Random(92)
    saw_path = saw_none = 0
    for _ in range(3000):
        n = rng.randint(1, 7)
        graph = {i: [] for i in range(n)}
        for _ in range(rng.randint(0, 10)):
            u = rng.randint(0, n - 1)
            v = rng.randint(0, n - 1)
            graph[u].append(v)
        trail = eulerian_path(graph, directed=True)
        total = sum(len(graph[u]) for u in graph)
        if has_eulerian_path(graph, True):
            assert trail is not None
            assert len(trail) == total + 1
            assert _edges_directed(graph) == _trail_edges_directed(trail)
            if has_eulerian_circuit(graph, True) and total > 0:
                assert trail[0] == trail[-1]
            saw_path += 1
        else:
            assert trail is None
            saw_none += 1
    assert saw_path > 0 and saw_none > 0


def test_undirected_square_circuit():
    g = {0: [1, 3], 1: [2], 2: [3]}  # edges 0-1, 0-3, 1-2, 2-3
    t = eulerian_path(g)
    assert t is not None
    assert t[0] == t[-1]
    assert _edges_undirected(g) == _trail_edges_undirected(t)
    assert has_eulerian_circuit(g, False)


def test_undirected_path_two_odd():
    g = {0: [1, 3], 1: [2], 2: [0]}  # degrees: 0->3, 1->2, 2->2, 3->1; odd = {0, 3}
    t = eulerian_path(g)
    assert t is not None
    assert t[0] in (0, 3) and t[-1] in (0, 3)
    assert _edges_undirected(g) == _trail_edges_undirected(t)
    assert has_eulerian_path(g, False)
    assert not has_eulerian_circuit(g, False)


def test_directed_circuit():
    g = {"a": ["b"], "b": ["c"], "c": ["a"]}
    t = eulerian_path(g, directed=True)
    assert t is not None
    assert t[0] == t[-1]
    assert _edges_directed(g) == _trail_edges_directed(t)


def test_directed_path():
    g = {"a": ["b"], "b": ["c"]}  # out-in diff: a +1, c -1 -> path a..c
    t = eulerian_path(g, directed=True)
    assert t == ["a", "b", "c"]


def test_star_three_odd_no_euler():
    g = {0: [1, 2, 3]}  # three odd-degree leaves + center degree 3 -> 4 odd vertices
    assert not has_eulerian_path(g, False)
    assert eulerian_path(g, False) is None


def test_disconnected_edges_no_euler():
    g = {0: [1], 2: [3]}  # two separate edges
    assert not has_eulerian_path(g, False)
    assert eulerian_path(g, False) is None


def test_empty_graph_is_circuit():
    assert has_eulerian_circuit({}, False)
    assert has_eulerian_circuit({}, True)


def test_single_edge_undirected():
    g = {0: [1]}
    t = eulerian_path(g)
    assert t == [0, 1] or t == [1, 0]
    assert has_eulerian_path(g, False)
    assert not has_eulerian_circuit(g, False)


def test_self_loop_undirected():
    g = {0: [0]}  # a single self-loop: degree 2, even -> circuit
    assert has_eulerian_circuit(g, False)
    t = eulerian_path(g)
    assert t == [0, 0]


def test_directed_unbalanced_no_euler():
    g = {"a": ["b", "c"], "b": [], "c": []}  # a out2 in0, others sink -> not Eulerian
    assert not has_eulerian_path(g, True)
    assert eulerian_path(g, True) is None


def test_start_vertex_respected_for_circuit():
    g = {0: [1, 3], 1: [2], 2: [3]}  # circuit
    t = eulerian_path(g, directed=False, start=2)
    assert t[0] == 2
    assert t[0] == t[-1]
    assert _edges_undirected(g) == _trail_edges_undirected(t)
