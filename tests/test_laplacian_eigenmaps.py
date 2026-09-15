import math
import random

import pytest

from quantforge import laplacian_eigenmaps, heat_knn_affinity


def _spearman(a, b):
    n = len(a)
    ra = sorted(range(n), key=lambda i: a[i])
    rb = sorted(range(n), key=lambda i: b[i])
    rka = [0] * n
    rkb = [0] * n
    for r, i in enumerate(ra):
        rka[i] = r
    for r, i in enumerate(rb):
        rkb[i] = r
    d2 = sum((rka[i] - rkb[i]) ** 2 for i in range(n))
    return 1 - 6 * d2 / (n * (n * n - 1))


def test_line_embedding_orders_by_position():
    X = [[i * 1.0, 0.3 * i] for i in range(25)]
    c = [cc[0] for cc in laplacian_eigenmaps(X, k=1, n_neighbors=4, t=5.0)["coords"]]
    assert abs(_spearman(c, list(range(25)))) > 0.95


def test_arc_embedding_orders_along_curve():
    th = [math.pi * i / 24 for i in range(25)]
    X = [[5 * math.cos(a), 5 * math.sin(a)] for a in th]
    c = [cc[0] for cc in laplacian_eigenmaps(X, k=1, n_neighbors=4, t=10.0)["coords"]]
    assert abs(_spearman(c, list(range(25)))) > 0.95


def test_two_clusters_separate():
    random.seed(0)
    A = [[random.gauss(0, 0.3), random.gauss(0, 0.3)] for _ in range(15)]
    B = [[random.gauss(10, 0.3), random.gauss(10, 0.3)] for _ in range(15)]
    c = [cc[0] for cc in laplacian_eigenmaps(A + B, k=2, n_neighbors=5, t=1.0)["coords"]]
    mid = (sum(c[:15]) / 15 + sum(c[15:]) / 15) / 2
    sA = set(1 if v > mid else 0 for v in c[:15])
    sB = set(1 if v > mid else 0 for v in c[15:])
    assert len(sA) == 1 and len(sB) == 1 and sA != sB


def test_affinity_symmetric_and_bounded():
    W = heat_knn_affinity([[0, 0], [1, 0], [2, 0]], 2, t=1.0)
    assert all(abs(W[i][j] - W[j][i]) < 1e-12 for i in range(3) for j in range(3))
    assert all(0 <= W[i][j] <= 1 for i in range(3) for j in range(3))


def test_errors():
    X = [[i, 0] for i in range(5)]
    with pytest.raises(ValueError):
        laplacian_eigenmaps(X, k=0)
    with pytest.raises(ValueError):
        laplacian_eigenmaps(X, k=5)     # k <= n-1 required
