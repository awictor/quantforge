"""Householder QR decomposition and QR least squares."""

import random

import pytest

from quantforge import qr_decomposition, qr_solve
from quantforge.ols import ols_fit


def test_reconstruction():
    rng = random.Random(1)
    m, n = 5, 3
    A = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(m)]
    Q, R = qr_decomposition(A)
    QR = [[sum(Q[i][k] * R[k][j] for k in range(m)) for j in range(n)]
          for i in range(m)]
    assert max(abs(QR[i][j] - A[i][j]) for i in range(m) for j in range(n)) < 1e-10


def test_q_orthogonal():
    rng = random.Random(2)
    m, n = 6, 4
    A = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(m)]
    Q, _ = qr_decomposition(A)
    QtQ = [[sum(Q[k][i] * Q[k][j] for k in range(m)) for j in range(m)]
           for i in range(m)]
    assert max(abs(QtQ[i][j] - (1.0 if i == j else 0.0))
               for i in range(m) for j in range(m)) < 1e-10


def test_r_upper_triangular():
    rng = random.Random(3)
    m, n = 5, 3
    A = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(m)]
    _, R = qr_decomposition(A)
    assert all(abs(R[i][j]) < 1e-12 for i in range(m) for j in range(min(i, n)))


def test_qr_solve_matches_ols():
    rng = random.Random(4)
    X = [[rng.gauss(0, 1), rng.gauss(0, 1)] for _ in range(50)]
    y = [2 * X[i][0] - 3 * X[i][1] + rng.gauss(0, 0.01) for i in range(50)]
    xq = qr_solve(X, y)
    xo = ols_fit(X, y, add_intercept=False)["coefficients"]
    assert all(abs(xq[k] - xo[k]) < 1e-8 for k in range(2))


def test_square_solve():
    A = [[2.0, 1.0], [1.0, 3.0]]
    b = [5.0, 10.0]
    x = qr_solve(A, b)
    # Exact solution of the 2x2 system.
    assert abs(A[0][0] * x[0] + A[0][1] * x[1] - b[0]) < 1e-10
    assert abs(A[1][0] * x[0] + A[1][1] * x[1] - b[1]) < 1e-10


def test_validation():
    with pytest.raises(ValueError):
        qr_decomposition([[1.0, 2.0]])          # m < n
    with pytest.raises(ValueError):
        qr_decomposition([])
    with pytest.raises(ValueError):
        qr_solve([[1.0], [2.0]], [1.0])         # b length mismatch
