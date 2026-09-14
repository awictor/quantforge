"""Tests for LCA (binary lifting), cross-checked against brute path-to-root ancestry."""

import random

import pytest

from quantforge.lca import LCA


def _build_random_tree(n, rng):
    adj = {0: []}
    for i in range(1, n):
        p = rng.randint(0, i - 1)
        adj.setdefault(i, [])
        adj.setdefault(p, [])
        adj[i].append(p)
        adj[p].append(i)
    return adj


def _brute_parent_depth(adj, root):
    par = {root: None}
    dep = {root: 0}
    st = [root]
    seen = {root}
    while st:
        u = st.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                par[w] = u
                dep[w] = dep[u] + 1
                st.append(w)
    return par, dep


def _brute_lca(par, u, v):
    au = set()
    x = u
    while x is not None:
        au.add(x)
        x = par[x]
    y = v
    while y not in au:
        y = par[y]
    return y


def test_fuzz_vs_brute():
    rng = random.Random(51)
    for _ in range(2000):
        n = rng.randint(1, 40)
        adj = _build_random_tree(n, rng)
        lca = LCA(adj, 0)
        par, dep = _brute_parent_depth(adj, 0)
        for _ in range(6):
            u = rng.randint(0, n - 1)
            v = rng.randint(0, n - 1)
            bw = _brute_lca(par, u, v)
            assert lca.query(u, v) == bw
            assert lca.depth(u) == dep[u]
            assert lca.distance(u, v) == dep[u] + dep[v] - 2 * dep[bw]


def _sample_tree():
    #        a
    #      /   \
    #     b     c
    #    / \     \
    #   d   e     f
    return {
        "a": ["b", "c"],
        "b": ["a", "d", "e"],
        "c": ["a", "f"],
        "d": ["b"],
        "e": ["b"],
        "f": ["c"],
    }


def test_sample_tree_lca():
    lca = LCA(_sample_tree(), "a")
    assert lca.query("d", "e") == "b"
    assert lca.query("d", "f") == "a"
    assert lca.query("d", "b") == "b"  # ancestor of itself
    assert lca.query("a", "f") == "a"


def test_sample_tree_depth():
    lca = LCA(_sample_tree(), "a")
    assert lca.depth("a") == 0
    assert lca.depth("b") == 1
    assert lca.depth("d") == 2


def test_sample_tree_distance():
    lca = LCA(_sample_tree(), "a")
    assert lca.distance("d", "f") == 4
    assert lca.distance("d", "e") == 2
    assert lca.distance("d", "d") == 0
    assert lca.distance("a", "d") == 2


def test_is_ancestor():
    lca = LCA(_sample_tree(), "a")
    assert lca.is_ancestor("a", "f") is True
    assert lca.is_ancestor("b", "d") is True
    assert lca.is_ancestor("b", "f") is False
    assert lca.is_ancestor("d", "d") is True  # self-ancestor


def test_query_is_symmetric():
    lca = LCA(_sample_tree(), "a")
    assert lca.query("d", "f") == lca.query("f", "d")
    assert lca.distance("d", "f") == lca.distance("f", "d")


def test_single_node_tree():
    lca = LCA({"x": []}, "x")
    assert lca.query("x", "x") == "x"
    assert lca.distance("x", "x") == 0
    assert lca.depth("x") == 0


def test_path_graph():
    # 0-1-2-3-4 line rooted at 0
    adj = {0: [1], 1: [0, 2], 2: [1, 3], 3: [2, 4], 4: [3]}
    lca = LCA(adj, 0)
    assert lca.query(2, 4) == 2
    assert lca.distance(0, 4) == 4
    assert lca.query(1, 3) == 1


def test_deep_tree_no_recursion_error():
    # a long chain would overflow recursive DFS; iterative build must handle it
    n = 5000
    adj = {0: []}
    for i in range(1, n):
        adj[i] = [i - 1]
        adj[i - 1].append(i)
    lca = LCA(adj, 0)
    assert lca.depth(n - 1) == n - 1
    assert lca.query(n - 1, n // 2) == n // 2
    assert lca.distance(0, n - 1) == n - 1


def test_bad_root_raises():
    with pytest.raises(ValueError):
        LCA({"a": []}, "z")


def test_disconnected_raises():
    with pytest.raises(ValueError):
        LCA({"a": ["b"], "b": ["a"], "c": []}, "a")


def test_root_from_middle():
    lca = LCA(_sample_tree(), "b")
    # rooting at b: a is a child of b, and f sits under c under a, so a is f's ancestor
    assert lca.query("d", "e") == "b"
    assert lca.depth("a") == 1
    assert lca.query("a", "f") == "a"
    assert lca.query("d", "f") == "b"  # d under b, f under a under b -> meet at b
