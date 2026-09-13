"""One-sided Jacobi SVD and the pseudo-inverse."""

import random

import pytest

from quantforge import svd, pseudo_inverse, qr_solve
from quantforge.pca import jacobi_eigen


def _matrix(seed=1, m=5, n=3):
    rng = random.Random(seed)
    return [[rng.gauss(0, 1) for _ in range(n)] for _ in range(m)]


def test_reconstruction():
    A = _matrix()
    m, n = len(A), len(A[0])
    U, s, V = svd(A)
    recon = [[sum(U[i][k] * s[k] * V[j][k] for k in range(n)) for j in range(n)]
             for i in range(m)]
    assert max(abs(recon[i][j] - A[i][j]) for i in range(m) for j in range(n)) < 1e-10


def test_singular_values_descending_nonneg():
    _, s, _ = svd(_matrix())
    assert all(s[i] >= s[i + 1] for i in range(len(s) - 1))
    assert all(x >= 0.0 for x in s)


def test_u_columns_orthonormal():
    A = _matrix()
    m, n = len(A), len(A[0])
    U, _, _ = svd(A)
    UtU = [[sum(U[k][i] * U[k][j] for k in range(m)) for j in range(n)]
           for i in range(n)]
    assert max(abs(UtU[i][j] - (1.0 if i == j else 0.0))
               for i in range(n) for j in range(n)) < 1e-10


def test_singular_values_match_eig_of_ata():
    A = _matrix()
    m, n = len(A), len(A[0])
    _, s, _ = svd(A)
    AtA = [[sum(A[k][i] * A[k][j] for k in range(m)) for j in range(n)]
           for i in range(n)]
    ev, _ = jacobi_eigen(AtA)
    for i in range(n):
        assert abs(s[i] * s[i] - ev[i]) < 1e-8


def test_pseudo_inverse_matches_qr_solve():
    A = _matrix()
    m, n = len(A), len(A[0])
    rng = random.Random(7)
    b = [rng.gauss(0, 1) for _ in range(m)]
    Ap = pseudo_inverse(A)
    xp = [sum(Ap[i][j] * b[j] for j in range(m)) for i in range(n)]
    xq = qr_solve(A, b)
    assert all(abs(xp[i] - xq[i]) < 1e-7 for i in range(n))


def test_validation():
    with pytest.raises(ValueError):
        svd([[1.0, 2.0]])          # m < n
    with pytest.raises(ValueError):
        svd([])
