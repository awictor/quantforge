import math
import random
from collections import defaultdict

import pytest

from quantforge import spectral_clustering, rbf_affinity
from quantforge.kmeans import kmeans


def _same_partition(labels, truth):
    g = defaultdict(set)
    t = defaultdict(set)
    for i, l in enumerate(labels):
        g[l].add(i)
    for i, x in enumerate(truth):
        t[x].add(i)
    return (sorted(map(frozenset, g.values()), key=min)
            == sorted(map(frozenset, t.values()), key=min))


def test_two_blobs():
    random.seed(0)
    A = [[random.gauss(0, 0.3), random.gauss(0, 0.3)] for _ in range(15)]
    B = [[random.gauss(5, 0.3), random.gauss(5, 0.3)] for _ in range(15)]
    r = spectral_clustering(A + B, 2, gamma=1.0, seed=1)
    assert _same_partition(r["labels"], [0] * 15 + [1] * 15)


def test_concentric_rings_where_kmeans_fails():
    random.seed(1)
    inner, outer = [], []
    for _ in range(30):
        a = random.uniform(0, 2 * math.pi)
        inner.append([math.cos(a) + random.gauss(0, 0.05), math.sin(a) + random.gauss(0, 0.05)])
        outer.append([4 * math.cos(a) + random.gauss(0, 0.05),
                      4 * math.sin(a) + random.gauss(0, 0.05)])
    X = inner + outer
    truth = [0] * 30 + [1] * 30
    assert _same_partition(spectral_clustering(X, 2, gamma=0.5, seed=2)["labels"], truth)
    assert not _same_partition(kmeans(X, 2, seed=2)["labels"], truth)


def test_three_blobs():
    random.seed(3)
    X, truth = [], []
    for c, (cx, cy) in enumerate([(0, 0), (8, 0), (4, 7)]):
        for _ in range(12):
            X.append([cx + random.gauss(0, 0.3), cy + random.gauss(0, 0.3)])
            truth.append(c)
    assert _same_partition(spectral_clustering(X, 3, gamma=0.5, seed=4)["labels"], truth)


def test_affinity_properties():
    W = rbf_affinity([[0, 0], [1, 0], [0, 1]], gamma=1.0)
    assert all(abs(W[i][j] - W[j][i]) < 1e-12 for i in range(3) for j in range(3))
    assert all(W[i][i] == 1.0 for i in range(3))
    assert all(0 < W[i][j] <= 1 for i in range(3) for j in range(3))


def test_invalid_k_raises():
    with pytest.raises(ValueError):
        spectral_clustering([[0, 0], [1, 1]], 0)
