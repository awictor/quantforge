"""Iterative linear solvers: CG, Gauss-Seidel, Jacobi."""

import random

import pytest

from quantforge import conjugate_gradient, gauss_seidel, jacobi
from quantforge.portopt import _invert


def _spd(n, seed):
    rng = random.Random(seed)
    A0 = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    return [[sum(A0[k][i] * A0[k][j] for k in range(n)) + (n if i == j else 0)
             for j in range(n)] for i in range(n)], [rng.gauss(0, 1) for _ in range(n)]


def test_cg_matches_direct_solve():
    n = 6
    A, b = _spd(n, 3)
    cg = conjugate_gradient(A, b)
    inv = _invert(A)
    xd = [sum(inv[i][j] * b[j] for j in range(n)) for i in range(n)]
    assert max(abs(cg["x"][i] - xd[i]) for i in range(n)) < 1e-8
    assert cg["n_iter"] <= n              # converges in <= n steps


def test_cg_identity():
    I = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    r = conjugate_gradient(I, [1, 2, 3, 4])
    assert r["x"] == [1.0, 2.0, 3.0, 4.0]
    assert r["n_iter"] == 1


def test_gauss_seidel_faster_than_jacobi():
    A = [[4, 1, 0], [1, 4, 1], [0, 1, 4]]
    b = [6, 12, 6]
    gs = gauss_seidel(A, b)
    jc = jacobi(A, b)
    inv = _invert(A)
    xd = [sum(inv[i][j] * b[j] for j in range(3)) for i in range(3)]
    assert max(abs(gs["x"][i] - xd[i]) for i in range(3)) < 1e-8
    assert max(abs(jc["x"][i] - xd[i]) for i in range(3)) < 1e-8
    assert gs["n_iter"] < jc["n_iter"]


def test_all_three_agree():
    A = [[4, 1, 0], [1, 4, 1], [0, 1, 4]]
    b = [6, 12, 6]
    cg = conjugate_gradient(A, b)["x"]
    gs = gauss_seidel(A, b)["x"]
    assert max(abs(cg[i] - gs[i]) for i in range(3)) < 1e-8


def test_large_spd_residual_small():
    A, b = _spd(30, 7)
    cg = conjugate_gradient(A, b)
    assert cg["residual_norm"] < 1e-8


def test_validation():
    with pytest.raises(ValueError):
        conjugate_gradient([[1, 2]], [1, 2])
    with pytest.raises(ValueError):
        gauss_seidel([[0, 1], [1, 0]], [1, 1])       # zero diagonal
