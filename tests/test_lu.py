"""LU decomposition, linear solve, and determinant."""

import random

import pytest

from quantforge import lu_decomposition, lu_solve, determinant


def test_known_determinants():
    assert abs(determinant([[1, 2], [3, 4]]) + 2.0) < 1e-12
    assert abs(determinant([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) - 1.0) < 1e-12
    assert abs(determinant([[2, 0, 0], [0, 3, 0], [0, 0, 4]]) - 24.0) < 1e-12
    assert abs(determinant([[6, 1, 1], [4, -2, 5], [2, 8, 7]]) + 306.0) < 1e-9


def test_singular_determinant_zero():
    assert determinant([[1, 2], [2, 4]]) == 0.0


def test_pa_equals_lu():
    rng = random.Random(1)
    n = 4
    A = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    L, U, piv, _ = lu_decomposition(A)
    PA = [A[piv[i]] for i in range(n)]
    LU = [[sum(L[i][k] * U[k][j] for k in range(n)) for j in range(n)]
          for i in range(n)]
    assert max(abs(PA[i][j] - LU[i][j]) for i in range(n) for j in range(n)) < 1e-10


def test_triangular_structure():
    rng = random.Random(2)
    n = 4
    A = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    L, U, _, _ = lu_decomposition(A)
    assert all(L[i][i] == 1.0 for i in range(n))
    assert all(abs(L[i][j]) < 1e-15 for i in range(n) for j in range(i + 1, n))
    assert all(abs(U[i][j]) < 1e-12 for i in range(n) for j in range(i))


def test_solve():
    rng = random.Random(3)
    n = 5
    A = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    b = [rng.gauss(0, 1) for _ in range(n)]
    x = lu_solve(A, b)
    resid = [sum(A[i][j] * x[j] for j in range(n)) - b[i] for i in range(n)]
    assert max(abs(r) for r in resid) < 1e-9


def test_pivoting_handles_zero_leading_pivot():
    # A[0][0] = 0 forces a row swap.
    A = [[0.0, 1.0], [1.0, 0.0]]
    x = lu_solve(A, [2.0, 3.0])
    assert abs(x[0] - 3.0) < 1e-12 and abs(x[1] - 2.0) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        lu_decomposition([[1.0, 2.0]])        # not square
    with pytest.raises(ValueError):
        lu_solve([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0])   # singular
    with pytest.raises(ValueError):
        lu_solve([[1.0, 0.0], [0.0, 1.0]], [1.0])        # b mismatch
