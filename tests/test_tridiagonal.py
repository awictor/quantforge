"""Tests for tridiagonal solvers, cross-checked against a dense Gaussian solve."""

import random

import pytest

from quantforge.tridiagonal import solve_tridiagonal, solve_cyclic_tridiagonal


def _dense_solve(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        M[col], M[piv] = M[piv], M[col]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col] / M[col][col]
                for k in range(col, n + 1):
                    M[r][k] -= f * M[col][k]
    return [M[i][n] / M[i][i] for i in range(n)]


def _close(a, b, tol=1e-6):
    return len(a) == len(b) and all(abs(x - y) <= tol * max(1, abs(x), abs(y)) for x, y in zip(a, b))


def test_fuzz_tridiagonal_vs_dense():
    rng = random.Random(341)
    for _ in range(3000):
        n = rng.randint(1, 20)
        diag = [rng.uniform(5, 10) * rng.choice([1, -1]) for _ in range(n)]
        lower = [0.0] + [rng.uniform(-2, 2) for _ in range(n - 1)]
        upper = [rng.uniform(-2, 2) for _ in range(n - 1)] + [0.0]
        rhs = [rng.uniform(-5, 5) for _ in range(n)]
        x = solve_tridiagonal(lower, diag, upper, rhs)
        A = [[0.0] * n for _ in range(n)]
        for i in range(n):
            A[i][i] = diag[i]
            if i > 0:
                A[i][i - 1] = lower[i]
            if i < n - 1:
                A[i][i + 1] = upper[i]
        assert _close(x, _dense_solve(A, rhs))


def test_fuzz_cyclic_vs_dense():
    rng = random.Random(342)
    for _ in range(3000):
        n = rng.randint(3, 20)
        diag = [rng.uniform(6, 12) * rng.choice([1, -1]) for _ in range(n)]
        lower = [rng.uniform(-2, 2) for _ in range(n)]
        upper = [rng.uniform(-2, 2) for _ in range(n)]
        rhs = [rng.uniform(-5, 5) for _ in range(n)]
        x = solve_cyclic_tridiagonal(lower, diag, upper, rhs)
        A = [[0.0] * n for _ in range(n)]
        for i in range(n):
            A[i][i] = diag[i]
            A[i][(i - 1) % n] += lower[i]
            A[i][(i + 1) % n] += upper[i]
        assert _close(x, _dense_solve(A, rhs), 1e-5)


def test_tridiagonal_explicit():
    # [[2,1],[1,2]] x = [3,3] -> [1,1]
    assert _close(solve_tridiagonal([0, 1], [2, 2], [1, 0], [3, 3]), [1, 1])


def test_identity():
    assert _close(solve_tridiagonal([0, 0, 0], [1, 1, 1], [0, 0, 0], [5, 6, 7]), [5, 6, 7])


def test_empty_and_single():
    assert solve_tridiagonal([], [], [], []) == []
    assert _close(solve_tridiagonal([0], [4], [0], [8]), [2])


def test_cyclic_reduces_to_known():
    # 3x3 circulant [[2,1,1],[1,2,1],[1,1,2]] x=[4,4,4] -> [1,1,1]
    x = solve_cyclic_tridiagonal([1, 1, 1], [2, 2, 2], [1, 1, 1], [4, 4, 4])
    assert _close(x, [1, 1, 1])


def test_zero_pivot_raises():
    with pytest.raises(ValueError):
        solve_tridiagonal([0], [0], [0], [1])


def test_cyclic_too_small_raises():
    with pytest.raises(ValueError):
        solve_cyclic_tridiagonal([1, 2], [1, 2], [1, 2], [1, 2])


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        solve_tridiagonal([0, 1], [1, 1], [1], [1, 1])


def test_solution_satisfies_system():
    # verify A x = rhs directly for a random case
    rng = random.Random(343)
    n = 8
    diag = [rng.uniform(5, 9) for _ in range(n)]
    lower = [0.0] + [rng.uniform(-1, 1) for _ in range(n - 1)]
    upper = [rng.uniform(-1, 1) for _ in range(n - 1)] + [0.0]
    rhs = [rng.uniform(-3, 3) for _ in range(n)]
    x = solve_tridiagonal(lower, diag, upper, rhs)
    for i in range(n):
        val = diag[i] * x[i]
        if i > 0:
            val += lower[i] * x[i - 1]
        if i < n - 1:
            val += upper[i] * x[i + 1]
        assert abs(val - rhs[i]) < 1e-9
