import random

import pytest

from quantforge import classical_mds, embedded_distances


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _maxdiff(A, B):
    n = len(A)
    return max(abs(A[i][j] - B[i][j]) for i in range(n) for j in range(n))


def test_recovers_2d_distances():
    random.seed(0)
    pts = [[random.uniform(-5, 5), random.uniform(-5, 5)] for _ in range(12)]
    D = embedded_distances(pts)
    res = classical_mds(D, 2)
    assert _maxdiff(D, embedded_distances(res["coords"])) < 1e-6


def test_rank_two_spectrum():
    random.seed(0)
    pts = [[random.uniform(-5, 5), random.uniform(-5, 5)] for _ in range(12)]
    res = classical_mds(embedded_distances(pts), 3)
    assert abs(res["eigenvalues"][2]) < 1e-6


def test_collinear_points_one_dimension():
    line = [[i * 1.5, 0.0] for i in range(8)]
    D = embedded_distances(line)
    res = classical_mds(D, 1)
    assert _maxdiff(D, embedded_distances(res["coords"])) < 1e-6


def test_equilateral_triangle():
    s = 1.0
    D = [[0, s, s], [s, 0, s], [s, s, 0]]
    res = classical_mds(D, 2)
    De = embedded_distances(res["coords"])
    for i in range(3):
        for j in range(3):
            if i != j:
                assert close(De[i][j], 1.0)


def test_recovers_3d_distances():
    random.seed(2)
    pts = [[random.uniform(-3, 3) for _ in range(3)] for _ in range(10)]
    D = embedded_distances(pts)
    res = classical_mds(D, 3)
    assert _maxdiff(D, embedded_distances(res["coords"])) < 1e-6


def test_errors():
    with pytest.raises(ValueError):
        classical_mds([[0, 1], [1, 0]], 0)
    with pytest.raises(ValueError):
        classical_mds([[0, 1], [1, 0]], 5)
