"""Power iteration and inverse iteration for eigenpairs."""

import random

import pytest

from quantforge import power_iteration, inverse_iteration, rayleigh_quotient
from quantforge.pca import jacobi_eigen


def _spd(n, seed):
    rng = random.Random(seed)
    A0 = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    return [[sum(A0[k][i] * A0[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]


def test_diagonal_dominant():
    r = power_iteration([[5, 0, 0], [0, 2, 0], [0, 0, 1]])
    assert abs(r["eigenvalue"] - 5.0) < 1e-8
    assert abs(r["eigenvector"][0] - 1.0) < 1e-6


def test_matches_jacobi_largest():
    A = _spd(5, 3)
    vals, _ = jacobi_eigen(A)
    r = power_iteration(A)
    assert abs(r["eigenvalue"] - max(vals)) < 1e-6


def test_eigenvector_residual():
    A = _spd(5, 3)
    r = power_iteration(A)
    v, lam = r["eigenvector"], r["eigenvalue"]
    Av = [sum(A[i][j] * v[j] for j in range(5)) for i in range(5)]
    assert max(abs(Av[i] - lam * v[i]) for i in range(5)) < 1e-5


def test_inverse_iteration_smallest():
    A = _spd(5, 3)
    vals, _ = jacobi_eigen(A)
    r = inverse_iteration(A, mu=0.0)
    assert abs(r["eigenvalue"] - min(vals)) < 1e-5


def test_inverse_iteration_nearest_shift():
    A = _spd(5, 3)
    vals, _ = jacobi_eigen(A)
    target = sorted(vals)[len(vals) // 2]
    r = inverse_iteration(A, mu=target + 0.01)
    assert abs(r["eigenvalue"] - target) < 1e-4


def test_rayleigh_exact_on_eigenvector():
    A = _spd(4, 5)
    vals, vecs = jacobi_eigen(A)
    assert abs(rayleigh_quotient(A, vecs[0]) - vals[0]) < 1e-8


def test_2x2_known():
    r = power_iteration([[2, 1], [1, 2]])       # eigenvalues 3, 1
    assert abs(r["eigenvalue"] - 3.0) < 1e-8


def test_validation():
    with pytest.raises(ValueError):
        rayleigh_quotient([[1, 0], [0, 1]], [0, 0])
