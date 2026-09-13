"""Non-negative least squares (Lawson-Hanson)."""

import random

import pytest

from quantforge import nnls


def test_exact_nonnegative_fit():
    A = [[1, 0], [0, 1], [1, 1]]
    b = [2, 3, 5]                         # x = [2, 3] fits exactly
    r = nnls(A, b)
    assert abs(r["x"][0] - 2) < 1e-6 and abs(r["x"][1] - 3) < 1e-6
    assert r["residual_norm"] < 1e-8


def test_negative_solution_clamped():
    A = [[1.0], [1.0]]
    b = [-2.0, -3.0]                      # OLS would give -2.5
    r = nnls(A, b)
    assert r["x"][0] == 0.0


def test_all_nonnegative():
    rng = random.Random(3)
    for _ in range(200):
        m, n = rng.randint(3, 8), rng.randint(1, 5)
        A = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(m)]
        b = [rng.gauss(0, 1) for _ in range(m)]
        assert all(xi > -1e-8 for xi in nnls(A, b)["x"])


def test_kkt_optimality():
    rng = random.Random(7)
    for _ in range(200):
        m, n = rng.randint(4, 10), rng.randint(2, 5)
        A = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(m)]
        b = [rng.gauss(0, 1) for _ in range(m)]
        x = nnls(A, b)["x"]
        resid = [b[i] - sum(A[i][j] * x[j] for j in range(n)) for i in range(m)]
        w = [sum(A[i][j] * resid[i] for i in range(m)) for j in range(n)]
        for j in range(n):
            if x[j] > 1e-6:
                assert abs(w[j]) < 1e-4        # passive: gradient ~ 0
            else:
                assert w[j] < 1e-4             # active: gradient <= 0


def test_beats_zero():
    A = [[1, 1], [1, 2], [1, 3]]
    b = [1, 2, 2]
    r = nnls(A, b)
    assert r["residual_norm"] < sum(bi * bi for bi in b) ** 0.5


def test_validation():
    with pytest.raises(ValueError):
        nnls([[1, 2]], [1, 2])            # row-count mismatch
    with pytest.raises(ValueError):
        nnls([], [])
