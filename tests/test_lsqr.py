import random

from quantforge import lsqr
from quantforge.svd import pseudo_inverse
from quantforge.lu import lu_solve


def close(a, b, tol=1e-5):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_overdetermined_matches_pseudo_inverse():
    random.seed(0)
    m, n = 15, 4
    A = [[random.gauss(0, 1) for _ in range(n)] for _ in range(m)]
    b = [random.gauss(0, 1) for _ in range(m)]
    res = lsqr(A, b, tol=1e-12)
    Ap = pseudo_inverse(A)
    xref = [sum(Ap[i][j] * b[j] for j in range(m)) for i in range(n)]
    for a, c in zip(res["x"], xref):
        assert close(a, c)


def test_square_matches_lu():
    random.seed(1)
    n = 8
    A = [[random.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        A[i][i] += n
    b = [random.gauss(0, 1) for _ in range(n)]
    res = lsqr(A, b, tol=1e-12)
    xref = lu_solve(A, b)
    for a, c in zip(res["x"], xref):
        assert close(a, c)


def test_consistent_system_zero_residual():
    random.seed(2)
    m, n = 10, 3
    A = [[random.gauss(0, 1) for _ in range(n)] for _ in range(m)]
    xtrue = [1.0, -2.0, 0.5]
    b = [sum(A[i][j] * xtrue[j] for j in range(n)) for i in range(m)]
    res = lsqr(A, b, tol=1e-12)
    assert res["residual"] < 1e-6
    for a, c in zip(res["x"], xtrue):
        assert close(a, c)


def test_polynomial_regression():
    random.seed(3)
    xs = [i * 0.4 for i in range(20)]
    ys = [2 - 1.5 * x + 0.3 * x * x + random.uniform(-0.1, 0.1) for x in xs]
    V = [[1.0, x, x * x] for x in xs]
    res = lsqr(V, ys, tol=1e-12)
    Ap = pseudo_inverse(V)
    ref = [sum(Ap[i][j] * ys[j] for j in range(20)) for i in range(3)]
    for a, c in zip(res["x"], ref):
        assert close(a, c, 1e-4)


def test_normal_equation_residual_is_small():
    random.seed(3)
    xs = [i * 0.4 for i in range(20)]
    ys = [2 - 1.5 * x + 0.3 * x * x for x in xs]
    V = [[1.0, x, x * x] for x in xs]
    res = lsqr(V, ys, tol=1e-12)
    r = [ys[i] - sum(V[i][j] * res["x"][j] for j in range(3)) for i in range(20)]
    atr = [sum(V[i][j] * r[i] for i in range(20)) for j in range(3)]
    assert all(abs(a) < 1e-4 for a in atr)
