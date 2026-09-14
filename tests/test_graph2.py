"""Graph algorithms II: union-find, MST, max flow."""

import itertools

import pytest

from quantforge import UnionFind, minimum_spanning_tree, max_flow


def test_union_find():
    uf = UnionFind()
    uf.union(1, 2)
    uf.union(3, 4)
    uf.union(2, 3)
    assert uf.connected(1, 4)
    assert not uf.connected(1, 5)
    assert uf.union(1, 4) is False       # already joined


def test_mst_weight_and_structure():
    nodes = ['A', 'B', 'C', 'D']
    edges = [('A', 'B', 1), ('A', 'C', 4), ('B', 'C', 2), ('B', 'D', 5), ('C', 'D', 1)]
    tree, total = minimum_spanning_tree(nodes, edges)
    assert total == 4
    assert len(tree) == len(nodes) - 1
    uf = UnionFind()
    for u, v, _ in tree:
        uf.union(u, v)
    assert len({uf.find(n) for n in nodes}) == 1     # spanning


def test_mst_matches_brute_force():
    nodes = ['A', 'B', 'C', 'D']
    edges = [('A', 'B', 1), ('A', 'C', 4), ('B', 'C', 2), ('B', 'D', 5), ('C', 'D', 1)]

    def brute(nodes, edges):
        best = None
        for combo in itertools.combinations(edges, len(nodes) - 1):
            uf = UnionFind()
            for x in nodes:
                uf.find(x)
            if all(uf.union(u, v) for u, v, _ in combo) and \
               len({uf.find(x) for x in nodes}) == 1:
                w = sum(w for _, _, w in combo)
                best = w if best is None else min(best, w)
        return best

    _, total = minimum_spanning_tree(nodes, edges)
    assert total == brute(nodes, edges)


def test_max_flow_clrs_classic():
    g = {'s': {'v1': 16, 'v2': 13}, 'v1': {'v2': 10, 'v3': 12},
         'v2': {'v1': 4, 'v4': 14}, 'v3': {'v2': 9, 't': 20},
         'v4': {'v3': 7, 't': 4}, 't': {}}
    assert max_flow(g, 's', 't') == 23


def test_max_flow_series_and_parallel():
    assert max_flow({'a': {'b': 5}, 'b': {'c': 3}, 'c': {}}, 'a', 'c') == 3
    assert max_flow({'s': {'a': 10, 'b': 10}, 'a': {'t': 10},
                     'b': {'t': 10}, 't': {}}, 's', 't') == 20


def test_validation():
    with pytest.raises(ValueError):
        max_flow({'a': {}}, 'a', 'a')                # source == sink
    with pytest.raises(ValueError):
        max_flow({'a': {'b': -1}, 'b': {}}, 'a', 'b')  # negative capacity
