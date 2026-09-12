"""K-means clustering."""

import random

import pytest

from quantforge import kmeans


def _blobs(seed=1):
    rng = random.Random(seed)
    X = []
    for cx, cy in [(0, 0), (10, 10), (0, 10)]:
        for _ in range(50):
            X.append([cx + rng.gauss(0, 0.5), cy + rng.gauss(0, 0.5)])
    return X


def test_separated_blobs_pure_clusters():
    X = _blobs()
    labels = kmeans(X, 3, seed=42)["labels"]
    # Each contiguous block of 50 came from one blob; it should get one label.
    assert all(len(set(labels[b * 50:(b + 1) * 50])) == 1 for b in range(3))


def test_centroids_near_true_centers():
    r = kmeans(_blobs(), 3, seed=42)
    cents = sorted(tuple(round(v) for v in c) for c in r["centroids"])
    assert cents == [(0, 0), (0, 10), (10, 10)]


def test_k_one_is_global_mean():
    X = _blobs()
    r = kmeans(X, 1)
    gm = [sum(p[j] for p in X) / len(X) for j in range(2)]
    assert all(abs(r["centroids"][0][j] - gm[j]) < 1e-9 for j in range(2))
    assert set(r["labels"]) == {0}


def test_reproducible_with_seed():
    X = _blobs()
    assert kmeans(X, 3, seed=7)["labels"] == kmeans(X, 3, seed=7)["labels"]


def test_more_clusters_lower_inertia():
    X = _blobs()
    assert kmeans(X, 5, seed=1)["inertia"] < kmeans(X, 2, seed=1)["inertia"]


def test_validation():
    X = _blobs()
    with pytest.raises(ValueError):
        kmeans(X, 0)
    with pytest.raises(ValueError):
        kmeans(X, len(X) + 1)
    with pytest.raises(ValueError):
        kmeans([], 3)
