import random

from quantforge import bicgstab, gmres
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
    A, b = _diag_dominant(10)
    xref = lu_solve(A, b)
    res = bicgstab(A, b, tol=1e-12)
    assert res["converged"]
    for a, c in zip(res["x"], xref):
        assert close(a, c)


def test_agrees_with_gmres():
    A, b = _diag_dominant(10)
    rb = bicgstab(A, b, tol=1e-12)
    rg = gmres(A, b, tol=1e-12)
    for a, c in zip(rb["x"], rg["x"]):
        assert close(a, c)


def test_matvec_callable():
    A, b = _diag_dominant(10)
    n = 10
    xref = lu_solve(A, b)
    mv = lambda v: [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]
    res = bicgstab(mv, b, tol=1e-12)
    for a, c in zip(res["x"], xref):
        assert close(a, c)


def test_spd_system():
    random.seed(2)
    n = 12
    M = [[random.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    S = [[sum(M[k][i] * M[k][j] for k in range(n)) + (n if i == j else 0) for j in range(n)]
         for i in range(n)]
    b = [random.gauss(0, 1) for _ in range(n)]
    xref = lu_solve(S, b)
    res = bicgstab(S, b, tol=1e-12)
    assert res["converged"]
    for a, c in zip(res["x"], xref):
        assert close(a, c, 1e-5)


def test_identity():
    I = [[1.0 if i == j else 0.0 for j in range(5)] for i in range(5)]
    res = bicgstab(I, [1, 2, 3, 4, 5])
    for a, c in zip(res["x"], [1, 2, 3, 4, 5]):
        assert close(a, c)


def test_large_diagonally_dominant():
    A, b = _diag_dominant(40, seed=3, boost=120)
    xref = lu_solve(A, b)
    res = bicgstab(A, b, tol=1e-10)
    assert res["converged"]
    assert max(abs(a - c) for a, c in zip(res["x"], xref)) < 1e-5
