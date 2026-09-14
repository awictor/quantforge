"""Tests for KDTree spatial queries, cross-checked against brute-force distances."""

import random

import pytest

from quantforge.kdtree import KDTree


def _bf_knn(pts, t, k):
    d = sorted((sum((a - b) ** 2 for a, b in zip(p, t)) ** 0.5, i) for i, p in enumerate(pts))
    return [(i, dist) for dist, i in d[:k]]


def _bf_radius(pts, t, r):
    out = [(i, sum((a - b) ** 2 for a, b in zip(p, t)) ** 0.5) for i, p in enumerate(pts)]
    return sorted([(i, dd) for i, dd in out if dd <= r + 1e-12], key=lambda x: (x[1], x[0]))


def _bf_range(pts, lo, hi):
    return sorted(
        i for i, p in enumerate(pts) if all(lo[a] <= p[a] <= hi[a] for a in range(len(lo)))
    )


def test_fuzz_knn_distance_multiset():
    rng = random.Random(31)
    for _ in range(1500):
        dim = rng.randint(1, 4)
        n = rng.randint(1, 40)
        pts = [tuple(rng.uniform(-10, 10) for _ in range(dim)) for _ in range(n)]
        tree = KDTree(pts)
        t = tuple(rng.uniform(-10, 10) for _ in range(dim))
        k = rng.randint(1, min(n, 6))
        got = tree.k_nearest(t, k)
        exp = _bf_knn(pts, t, k)
        # distances must match as a sorted multiset (index ties may resolve differently)
        assert sorted(round(d, 9) for _, d in got) == sorted(round(d, 9) for _, d in exp)


def test_fuzz_nearest_distance():
    rng = random.Random(32)
    for _ in range(1000):
        dim = rng.randint(1, 3)
        n = rng.randint(1, 30)
        pts = [tuple(rng.uniform(-5, 5) for _ in range(dim)) for _ in range(n)]
        tree = KDTree(pts)
        t = tuple(rng.uniform(-5, 5) for _ in range(dim))
        assert tree.nearest(t)[1] == pytest.approx(_bf_knn(pts, t, 1)[0][1], abs=1e-12)


def test_fuzz_within_radius():
    rng = random.Random(33)
    for _ in range(1500):
        dim = rng.randint(1, 4)
        n = rng.randint(1, 40)
        pts = [tuple(rng.uniform(-10, 10) for _ in range(dim)) for _ in range(n)]
        tree = KDTree(pts)
        t = tuple(rng.uniform(-10, 10) for _ in range(dim))
        r = rng.uniform(0, 8)
        got = tree.within_radius(t, r)
        exp = _bf_radius(pts, t, r)
        assert sorted(i for i, _ in got) == sorted(i for i, _ in exp)


def test_fuzz_range_search():
    rng = random.Random(34)
    for _ in range(1500):
        dim = rng.randint(1, 4)
        n = rng.randint(1, 40)
        pts = [tuple(rng.uniform(-10, 10) for _ in range(dim)) for _ in range(n)]
        tree = KDTree(pts)
        lo = tuple(rng.uniform(-10, 0) for _ in range(dim))
        hi = tuple(l + rng.uniform(0, 15) for l in lo)
        assert tree.range_search(lo, hi) == _bf_range(pts, lo, hi)


def test_knn_returns_sorted_by_distance():
    pts = [(0, 0), (1, 0), (0, 1), (5, 5), (2, 2)]
    tree = KDTree(pts)
    res = tree.k_nearest((0, 0), 5)
    dists = [d for _, d in res]
    assert dists == sorted(dists)


def test_index_maps_to_distance():
    pts = [(0, 0), (3, 4), (1, 1)]
    tree = KDTree(pts)
    idx, dist = tree.nearest((0, 0))
    assert idx == 0
    assert dist == pytest.approx(0.0)
    idx, dist = tree.nearest((3, 4))
    assert idx == 1
    assert dist == pytest.approx(0.0)


def test_within_radius_inclusive():
    pts = [(0, 0), (3, 0), (5, 0)]
    tree = KDTree(pts)
    res = tree.within_radius((0, 0), 3.0)
    assert sorted(i for i, _ in res) == [0, 1]  # distance exactly 3 is included


def test_range_search_inclusive_box():
    pts = [(0, 0), (1, 1), (2, 2), (3, 3)]
    tree = KDTree(pts)
    assert tree.range_search((1, 1), (2, 2)) == [1, 2]


def test_single_point():
    tree = KDTree([(2.0, 3.0)])
    assert tree.nearest((0, 0))[0] == 0
    assert tree.k_nearest((0, 0), 5) == [(0, pytest.approx(13 ** 0.5))]


def test_one_dimensional():
    pts = [(x,) for x in (5, 1, 9, 3, 7)]
    tree = KDTree(pts)
    idx, dist = tree.nearest((4,))
    assert pts[idx] in {(3,), (5,)}  # both are distance 1 from 4 (a tie)
    assert dist == pytest.approx(1.0)
    # unambiguous target
    idx2, dist2 = tree.nearest((8.6,))
    assert pts[idx2] == (9,)
    assert dist2 == pytest.approx(0.4)


def test_len():
    assert len(KDTree([(1, 2), (3, 4), (5, 6)])) == 3
    assert len(KDTree([])) == 0


def test_k_larger_than_n_returns_all():
    pts = [(0, 0), (1, 1)]
    tree = KDTree(pts)
    res = tree.k_nearest((0, 0), 10)
    assert len(res) == 2


def test_empty_tree_queries():
    tree = KDTree([])
    assert tree.k_nearest((0, 0), 3) == []
    assert tree.within_radius((0, 0), 5) == []


def test_k_zero_returns_empty():
    tree = KDTree([(0, 0), (1, 1)])
    assert tree.k_nearest((0, 0), 0) == []


def test_empty_nearest_raises():
    with pytest.raises(ValueError):
        KDTree([]).nearest((0, 0))


def test_negative_k_raises():
    with pytest.raises(ValueError):
        KDTree([(0, 0)]).k_nearest((0, 0), -1)


def test_negative_radius_raises():
    with pytest.raises(ValueError):
        KDTree([(0, 0)]).within_radius((0, 0), -1)


def test_dimension_mismatch_raises():
    with pytest.raises(ValueError):
        KDTree([(1, 2), (3,)])


def test_range_bounds_dimension_mismatch_raises():
    tree = KDTree([(1, 2), (3, 4)])
    with pytest.raises(ValueError):
        tree.range_search((0,), (1, 1))
