"""Tests for DAG longest/shortest paths and reachability, cross-checked by path enumeration."""

import random

import pytest

from quantforge.dag_paths import dag_shortest_path, dag_longest_path, transitive_closure


def _all_paths(graph, s, t):
    res = []

    def dfs(u, acc, seen):
        if u == t:
            res.append(acc)
            return
        for v, w in graph[u]:
            if v not in seen:
                dfs(v, acc + w, seen | {v})

    dfs(s, 0.0, {s})
    return res


def _random_dag(rng, n):
    graph = {i: [] for i in range(n)}
    for u in range(n):
        for v in range(u + 1, n):  # only forward edges -> acyclic
            if rng.random() < 0.4:
                graph[u].append((v, rng.randint(-5, 5)))
    return graph


def test_fuzz_short_long_vs_enumeration():
    rng = random.Random(401)
    for _ in range(3000):
        n = rng.randint(1, 8)
        graph = _random_dag(rng, n)
        s = rng.randint(0, n - 1)
        dshort = dag_shortest_path(graph, s)
        dlong = dag_longest_path(graph, s)
        for t in range(n):
            paths = _all_paths(graph, s, t)
            if t == s:
                assert dshort[t] == min([0.0] + paths)
                assert dlong[t] == max([0.0] + paths)
            elif paths:
                assert dshort[t] == min(paths)
                assert dlong[t] == max(paths)
            else:
                assert dshort[t] == float("inf")
                assert dlong[t] == float("-inf")


def test_fuzz_path_reconstruction():
    rng = random.Random(402)
    for _ in range(1000):
        n = rng.randint(1, 8)
        graph = _random_dag(rng, n)
        s = rng.randint(0, n - 1)
        for t in range(n):
            d, p = dag_shortest_path(graph, s, t)
            if d != float("inf"):
                assert p[0] == s and p[-1] == t
                w = sum(
                    dict(graph[p[i]])[p[i + 1]] for i in range(len(p) - 1)
                )
                assert abs(w - d) < 1e-9


def test_fuzz_transitive_closure():
    rng = random.Random(403)
    for _ in range(2000):
        n = rng.randint(1, 8)
        graph = {i: [] for i in range(n)}
        for u in range(n):
            for v in range(n):
                if u != v and rng.random() < 0.3:
                    graph[u].append(v)
        tc = transitive_closure(graph)
        for s in range(n):
            seen = set()
            st = [s]
            while st:
                u = st.pop()
                for v in graph[u]:
                    if v not in seen:
                        seen.add(v)
                        st.append(v)
            assert tc[s] == seen


def test_critical_path_explicit():
    g = {0: [(1, 3), (2, 2)], 1: [(3, 4)], 2: [(3, 1)], 3: [(4, 2)], 4: []}
    d, p = dag_longest_path(g, 0, 4)
    assert d == 9
    assert p == [0, 1, 3, 4]


def test_shortest_path_explicit():
    g = {0: [(1, 3), (2, 2)], 1: [(3, 4)], 2: [(3, 1)], 3: [(4, 2)], 4: []}
    d, p = dag_shortest_path(g, 0, 4)
    assert d == 5
    assert p == [0, 2, 3, 4]


def test_negative_weights():
    g = {0: [(1, 5), (2, 2)], 1: [(3, -10)], 2: [(3, 1)], 3: []}
    assert dag_shortest_path(g, 0)[3] == -5  # 0->1->3 = 5 - 10
    assert dag_longest_path(g, 0)[3] == 3     # 0->2->3 = 2 + 1


def test_unreachable():
    g = {0: [(1, 1)], 1: [], 2: []}
    assert dag_shortest_path(g, 0)[2] == float("inf")
    assert dag_longest_path(g, 0)[2] == float("-inf")
    assert dag_shortest_path(g, 0, 2) == (float("inf"), [])


def test_transitive_closure_plain_adjacency():
    g = {"a": ["b"], "b": ["c"], "c": []}
    tc = transitive_closure(g)
    assert tc["a"] == {"b", "c"}
    assert tc["b"] == {"c"}
    assert tc["c"] == set()


def test_transitive_closure_with_cycle():
    g = {0: [1], 1: [2], 2: [0]}  # cycle allowed for closure
    tc = transitive_closure(g)
    assert tc[0] == {0, 1, 2}


def test_single_node():
    assert dag_shortest_path({0: []}, 0) == {0: 0.0}
    assert dag_longest_path({0: []}, 0) == {0: 0.0}


def test_cycle_raises():
    with pytest.raises(ValueError):
        dag_shortest_path({0: [(1, 1)], 1: [(0, 1)]}, 0)


def test_source_not_in_graph_raises():
    with pytest.raises(ValueError):
        dag_shortest_path({0: []}, 5)
