"""Exact integer/rational linear algebra (Bareiss fraction-free elimination)."""

import random
from fractions import Fraction

import pytest

from quantforge import bareiss_determinant, rational_solve, rational_inverse


def _cofactor_det(m):
    n = len(m)
    if n == 1:
        return m[0][0]
    s = 0
    for j in range(n):
        minor = [[m[i][k] for k in range(n) if k != j] for i in range(1, n)]
        s += ((-1) ** j) * m[0][j] * _cofactor_det(minor)
    return s


def test_determinant_matches_cofactor():
    rng = random.Random(1)
    for _ in range(200):
        n = rng.randint(1, 5)
        m = [[rng.randint(-9, 9) for _ in range(n)] for _ in range(n)]
        assert bareiss_determinant(m) == _cofactor_det(m)


def test_determinant_known():
    assert bareiss_determinant([[1, 2], [3, 4]]) == -2
    assert bareiss_determinant([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) == 1
    assert bareiss_determinant([[1, 2], [2, 4]]) == 0        # singular


def test_rational_solve_exact():
    rng = random.Random(2)
    solved = 0
    for _ in range(200):
        n = rng.randint(1, 5)
        A = [[rng.randint(-9, 9) for _ in range(n)] for _ in range(n)]
        if bareiss_determinant(A) == 0:
            continue
        b = [rng.randint(-9, 9) for _ in range(n)]
        x = rational_solve(A, b)
        for i in range(n):
            assert sum(A[i][j] * x[j] for j in range(n)) == Fraction(b[i])
        solved += 1
    assert solved > 100


def test_hilbert_matrix_solve_is_exact():
    n = 6
    H = [[Fraction(1, i + j + 1) for j in range(n)] for i in range(n)]
    b = [sum(H[i]) for i in range(n)]        # true solution is all ones
    x = rational_solve(H, b)
    assert all(v == 1 for v in x)


def test_inverse_times_matrix_is_identity():
    A = [[2, 1, 1], [1, 3, 2], [1, 0, 0]]
    inv = rational_inverse(A)
    n = len(A)
    for i in range(n):
        for j in range(n):
            s = sum(inv[i][k] * A[k][j] for k in range(n))
            assert s == (1 if i == j else 0)


def test_singular_raises():
    with pytest.raises(ValueError):
        rational_solve([[1, 2], [2, 4]], [1, 2])
    with pytest.raises(ValueError):
        rational_inverse([[1, 2], [2, 4]])


def test_validation():
    with pytest.raises(ValueError):
        bareiss_determinant([])
    with pytest.raises(ValueError):
        bareiss_determinant([[1, 2, 3], [4, 5, 6]])      # not square
    with pytest.raises(ValueError):
        rational_solve([[1, 2], [3, 4]], [1])            # b length mismatch
