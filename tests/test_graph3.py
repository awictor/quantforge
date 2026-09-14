"""Graph centrality: PageRank, degree, closeness, betweenness."""

import pytest

from quantforge import (
    pagerank,
    degree_centrality,
    closeness_centrality,
    betweenness_centrality,
)

STAR = {0: [1, 2, 3, 4], 1: [0], 2: [0], 3: [0], 4: [0]}
PATH = {0: [1], 1: [0, 2], 2: [1, 3], 3: [2, 4], 4: [3]}


def test_pagerank_sums_to_one():
    g = {'A': ['B', 'C'], 'B': ['C'], 'C': ['A']}
    pr = pagerank(g)
    assert abs(sum(pr.values()) - 1.0) < 1e-9


def test_pagerank_symmetric_ring_is_uniform():
    ring = {0: [1], 1: [2], 2: [3], 3: [0]}
    pr = pagerank(ring)
    assert max(pr.values()) - min(pr.values()) < 1e-6


def test_pagerank_handles_dangling_node():
    pr = pagerank({'A': ['B'], 'B': ['C'], 'C': []})
    assert abs(sum(pr.values()) - 1.0) < 1e-7


def test_degree_centrality_star():
    dc = degree_centrality(STAR)
    assert dc[0] == 1.0
    assert all(abs(dc[i] - 0.25) < 1e-12 for i in (1, 2, 3, 4))


def test_closeness_centrality_star():
    cc = closeness_centrality(STAR)
    assert cc[0] == max(cc.values())            # center is closest to all


def test_betweenness_star_and_path():
    bs = betweenness_centrality(STAR)
    assert bs[0] == max(bs.values()) and bs[0] > 0    # center on every leaf-leaf path
    assert all(bs[i] == 0.0 for i in (1, 2, 3, 4))
    bp = betweenness_centrality(PATH)
    assert bp[2] == max(bp.values())            # middle of the path


def test_validation():
    with pytest.raises(ValueError):
        pagerank({})
    with pytest.raises(ValueError):
        pagerank({'A': []}, damping=1.5)
