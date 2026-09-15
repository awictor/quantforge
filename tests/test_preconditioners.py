import random

from quantforge import (
    jacobi_preconditioner,
    incomplete_cholesky,
    ic_apply,
    preconditioned_cg,
)
from quantforge.lu import lu_solve


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _spd(n, seed=0, jitter=None):
    random.seed(seed)
    M = [[random.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    A = [[sum(M[k][i] * M[k][j] for k in range(n)) + (n if i == j else 0.0) for j in range(n)]
         for i in range(n)]
    if jitter:
        for i in range(n):
            A[i][i] += jitter(i)
    b = [random.gauss(0, 1) for _ in range(n)]
    return A, b


def test_plain_pcg_matches_lu():
    A, b = _spd(12)
    xref = lu_solve(A, b)
    r = preconditioned_cg(A, b, tol=1e-12)
    assert r["converged"]
    for a, c in zip(r["x"], xref):
        assert close(a, c)


def test_jacobi_preconditioner():
    A, b = _spd(12)
    xref = lu_solve(A, b)
    r = preconditioned_cg(A, b, apply_minv=jacobi_preconditioner(A), tol=1e-12)
    assert r["converged"]
    for a, c in zip(r["x"], xref):
        assert close(a, c)


def test_ic0_factor_matches_dense_cholesky():
    A, _ = _spd(12)
    n = len(A)
    L = incomplete_cholesky(A)
    LLt = [[sum(L[i][k] * L[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
    assert max(abs(A[i][j] - LLt[i][j]) for i in range(n) for j in range(n)) < 1e-8


def test_ic_preconditioned_cg():
    A, b = _spd(12)
    xref = lu_solve(A, b)
    L = incomplete_cholesky(A)
    r = preconditioned_cg(A, b, apply_minv=lambda v: ic_apply(L, v), tol=1e-12)
    assert r["converged"]
    for a, c in zip(r["x"], xref):
        assert close(a, c)


def test_ic_reduces_iterations_on_ill_conditioned():
    A, b = _spd(30, seed=2, jitter=lambda i: 10 ** (i % 5))
    plain = preconditioned_cg(A, b, tol=1e-8)
    L = incomplete_cholesky(A)
    icp = preconditioned_cg(A, b, apply_minv=lambda v: ic_apply(L, v), tol=1e-8)
    assert plain["converged"] and icp["converged"]
    assert icp["n_iter"] <= plain["n_iter"]


def test_ic0_preserves_sparsity():
    A = [[0.0] * 5 for _ in range(5)]
    for i in range(5):
        A[i][i] = 4.0
    for i in range(4):
        A[i][i + 1] = A[i + 1][i] = -1.0     # tridiagonal
    L = incomplete_cholesky(A)
    for i in range(5):
        for j in range(i - 1):               # below the first sub-diagonal must stay zero
            assert L[i][j] == 0.0
