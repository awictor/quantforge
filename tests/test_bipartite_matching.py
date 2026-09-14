"""Tests for Hopcroft-Karp bipartite matching, cross-checked against brute Kuhn matching."""

import random

from quantforge.bipartite_matching import (
    maximum_bipartite_matching,
    maximum_matching_size,
    minimum_vertex_cover,
)


def _brute_matching_size(adj):
    mr = {}

    def try_k(u, seen):
        for v in adj[u]:
            if v in seen:
                continue
            seen.add(v)
            if v not in mr or try_k(mr[v], seen):
                mr[v] = u
                return True
        return False

    cnt = 0
    for u in adj:
        if try_k(u, set()):
            cnt += 1
    return cnt


def _valid_matching(adj, m):
    rights = list(m.values())
    if len(rights) != len(set(rights)):
        return False
    return all(v in adj[u] for u, v in m.items())


def _valid_cover(adj, lc, rc):
    return all(u in lc or v in rc for u in adj for v in adj[u])


def _random_graph(rng):
    nl = rng.randint(1, 8)
    nr = rng.randint(1, 8)
    adj = {}
    for u in range(nl):
        adj[f"L{u}"] = [f"R{v}" for v in range(nr) if rng.random() < 0.35]
    return adj


def test_fuzz_matching_size_vs_brute():
    rng = random.Random(81)
    for _ in range(4000):
        adj = _random_graph(rng)
        m = maximum_bipartite_matching(adj)
        assert _valid_matching(adj, m)
        assert len(m) == _brute_matching_size(adj)


def test_fuzz_konig_identity():
    rng = random.Random(82)
    for _ in range(4000):
        adj = _random_graph(rng)
        m = maximum_bipartite_matching(adj)
        lc, rc = minimum_vertex_cover(adj)
        assert len(lc) + len(rc) == len(m)  # Konig's theorem
        assert _valid_cover(adj, lc, rc)


def test_explicit_perfect_matching():
    adj = {"a": ["x", "y"], "b": ["x"], "c": ["y", "z"]}
    m = maximum_bipartite_matching(adj)
    assert len(m) == 3
    assert _valid_matching(adj, m)


def test_k22_complete():
    adj = {"a": ["x", "y"], "b": ["x", "y"]}
    assert maximum_matching_size(adj) == 2


def test_star_left():
    assert maximum_matching_size({"a": ["x", "y", "z"]}) == 1


def test_star_right():
    assert maximum_matching_size({"a": ["x"], "b": ["x"], "c": ["x"]}) == 1


def test_cover_matches_matching_explicit():
    adj = {"a": ["x", "y"], "b": ["x"], "c": ["y", "z"]}
    lc, rc = minimum_vertex_cover(adj)
    assert len(lc) + len(rc) == 3
    assert _valid_cover(adj, lc, rc)


def test_empty_graph():
    assert maximum_bipartite_matching({}) == {}
    assert maximum_matching_size({}) == 0
    assert minimum_vertex_cover({}) == (set(), set())


def test_isolated_left_vertex():
    assert maximum_matching_size({"a": []}) == 0
    assert minimum_vertex_cover({"a": []}) == (set(), set())


def test_overlapping_left_right_labels():
    # left and right share the value "1" but live in separate namespaces
    adj = {1: [1, 2], 2: [1]}
    m = maximum_bipartite_matching(adj)
    assert len(m) == 2
    assert _valid_matching(adj, m)


def test_matching_is_dict_left_to_right():
    adj = {"a": ["x"], "b": ["y"]}
    m = maximum_bipartite_matching(adj)
    assert m == {"a": "x", "b": "y"}


def test_disjoint_components():
    adj = {"a": ["x"], "b": ["y"], "c": ["z"]}
    assert maximum_matching_size(adj) == 3
