import math
import random

import pytest

from quantforge import isomap, knn_graph
from quantforge.mds import embedded_distances


def test_flat_plane_preserves_distances():
    random.seed(0)
    pts2d = [[random.uniform(0, 5), random.uniform(0, 5)] for _ in range(40)]
    X = [[p[0], p[1], 0.5 * p[0] + 0.3 * p[1]] for p in pts2d]  # a plane in 3D
    res = isomap(X, k=2, n_neighbors=10)
    D3 = embedded_distances(X)
    De = embedded_distances(res["coords"])
    num = sum(abs(D3[i][j] - De[i][j]) for i in range(40) for j in range(i + 1, 40))
    den = sum(D3[i][j] for i in range(40) for j in range(i + 1, 40))
    assert num / den < 0.05


def test_arc_unrolls_monotonically():
    th = [math.pi * i / 29 for i in range(30)]
    X = [[5 * math.cos(a), 5 * math.sin(a)] for a in th]
    res = isomap(X, k=1, n_neighbors=3)
    c = [cc[0] for cc in res["coords"]]
    order = sorted(range(30), key=lambda i: c[i])
    assert order == list(range(30)) or order[::-1] == list(range(30))


def test_geodesic_exceeds_euclidean():
    th = [math.pi * i / 29 for i in range(30)]
    X = [[5 * math.cos(a), 5 * math.sin(a)] for a in th]
    res = isomap(X, k=1, n_neighbors=3)
    c = [cc[0] for cc in res["coords"]]
    euc = math.sqrt(sum((X[0][d] - X[29][d]) ** 2 for d in range(2)))
    assert abs(c[0] - c[29]) > euc      # half-circle arc length pi*r ~ 15.7 > chord 10


def test_knn_graph_symmetric():
    g = knn_graph([[0, 0], [1, 0], [2, 0], [3, 0]], 2)
    for i in g:
        for j in g[i]:
            assert i in g[j]
            assert abs(g[i][j] - g[j][i]) < 1e-12


def test_disconnected_raises():
    with pytest.raises(ValueError):
        isomap([[0, 0], [0, 1], [100, 100], [100, 101]], k=1, n_neighbors=1)
