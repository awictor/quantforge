"""Agglomerative hierarchical clustering."""

import random

import pytest

from quantforge import linkage, fcluster


def _blobs(seed=1):
    rng = random.Random(seed)
    X = []
    for cx, cy in [(0, 0), (10, 10), (0, 10)]:
        for _ in range(15):
            X.append([cx + rng.gauss(0, 0.3), cy + rng.gauss(0, 0.3)])
    return X


def test_merge_count_and_monotone():
    X = _blobs()
    mg = linkage(X, "average")
    assert len(mg) == len(X) - 1
    ds = [m[2] for m in mg]
    assert all(ds[i] <= ds[i + 1] + 1e-9 for i in range(len(ds) - 1))


def test_cut_into_three_pure_clusters():
    X = _blobs()
    labels = fcluster(X, linkage(X, "average"), n_clusters=3)
    assert len(set(labels)) == 3
    assert all(len(set(labels[b * 15:(b + 1) * 15])) == 1 for b in range(3))


def test_distance_threshold_cut():
    X = _blobs()
    labels = fcluster(X, linkage(X, "average"), distance_threshold=5.0)
    assert len(set(labels)) == 3


def test_extreme_cuts():
    X = _blobs()
    mg = linkage(X, "average")
    assert len(set(fcluster(X, mg, n_clusters=1))) == 1
    assert len(set(fcluster(X, mg, n_clusters=len(X)))) == len(X)


def test_complete_at_least_single_final_merge():
    X = _blobs()
    assert linkage(X, "complete")[-1][2] >= linkage(X, "single")[-1][2]


def test_validation():
    X = _blobs()
    mg = linkage(X, "average")
    with pytest.raises(ValueError):
        linkage([[1.0]])                      # < 2 points
    with pytest.raises(ValueError):
        linkage(X, "bad")                     # unknown method
    with pytest.raises(ValueError):
        fcluster(X, mg)                       # no cut criterion
    with pytest.raises(ValueError):
        fcluster(X, mg, n_clusters=1, distance_threshold=1.0)  # both given
