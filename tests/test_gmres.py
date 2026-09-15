import random

from quantforge import gmres
from quantforge.lu import lu_solve


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _diag_dominant(n, seed=0, boost=None):
    random.seed(seed)
    A = [[random.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        A[i][i] += (boost if boost is not None else n)
    b = [random.gauss(0, 1) for _ in range(n)]
    return A, b


def test_nonsymmetric_solve_matches_lu():
    A, b = _diag_dominant(8)
    xref = lu_solve(A, b)
    res = gmres(A, b, tol=1e-12)
    assert res["converged"]
    for a, c in zip(res["x"], xref):
        assert close(a, c)


def test_residual_monotone():
    A, b = _diag_dominant(8)
    res = gmres(A, b, tol=1e-12)
    rs = res["residuals"]
    assert all(rs[i + 1] <= rs[i] + 1e-12 for i in range(len(rs) - 1))


def test_matvec_callable():
    A, b = _diag_dominant(8)
    n = 8
    xref = lu_solve(A, b)
    mv = lambda v: [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]
    res = gmres(mv, b, tol=1e-12)
    for a, c in zip(res["x"], xref):
        assert close(a, c)


def test_finite_termination():
    A, b = _diag_dominant(8)
    res = gmres(A, b, tol=1e-14)
    assert res["n_iter"] <= 8


def test_restart():
    A, b = _diag_dominant(20, seed=2, boost=40)
    xref = lu_solve(A, b)
    res = gmres(A, b, tol=1e-10, restart=10)
    assert res["converged"]
    for a, c in zip(res["x"], xref):
        assert close(a, c, 1e-5)


def test_identity():
    I = [[1.0 if i == j else 0.0 for j in range(5)] for i in range(5)]
    res = gmres(I, [1, 2, 3, 4, 5])
    for a, c in zip(res["x"], [1, 2, 3, 4, 5]):
        assert close(a, c)
