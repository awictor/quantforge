"""Matrix exponential (scaling-and-squaring with Pade)."""

import math
import random

import pytest

from quantforge import matrix_exp


def test_exp_zero_is_identity():
    r = matrix_exp([[0.0, 0.0], [0.0, 0.0]])
    assert r == [[1.0, 0.0], [0.0, 1.0]]


def test_exp_diagonal():
    r = matrix_exp([[1.0, 0.0], [0.0, 2.0]])
    assert abs(r[0][0] - math.e) < 1e-10
    assert abs(r[1][1] - math.e ** 2) < 1e-10
    assert abs(r[0][1]) < 1e-12


def test_exp_nilpotent():
    # N = [[0,1],[0,0]] is nilpotent: exp(N) = I + N.
    r = matrix_exp([[0.0, 1.0], [0.0, 0.0]])
    assert abs(r[0][0] - 1.0) < 1e-12
    assert abs(r[0][1] - 1.0) < 1e-12
    assert abs(r[1][1] - 1.0) < 1e-12


def test_inverse_relation():
    rng = random.Random(1)
    A = [[rng.gauss(0, 1) for _ in range(3)] for _ in range(3)]
    negA = [[-A[i][j] for j in range(3)] for i in range(3)]
    eA = matrix_exp(A)
    emA = matrix_exp(negA)
    prod = [[sum(eA[i][k] * emA[k][j] for k in range(3)) for j in range(3)]
            for i in range(3)]
    assert max(abs(prod[i][j] - (1.0 if i == j else 0.0))
               for i in range(3) for j in range(3)) < 1e-9


def test_markov_generator_gives_stochastic_matrix():
    # Rows of a generator sum to zero -> exp rows sum to one, entries non-negative.
    Q = [[-0.3, 0.2, 0.1], [0.1, -0.2, 0.1], [0.05, 0.05, -0.1]]
    P = matrix_exp(Q)
    for row in P:
        assert abs(sum(row) - 1.0) < 1e-9
    assert all(v >= -1e-12 for row in P for v in row)


def test_matches_series_for_small_norm():
    A = [[0.1, 0.05], [-0.02, 0.08]]
    # exp(A) ~ I + A + A^2/2 + A^3/6 for small A.
    n = 2
    I = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    def mm(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)]
                for i in range(n)]
    A2 = mm(A, A)
    A3 = mm(A2, A)
    approx = [[I[i][j] + A[i][j] + A2[i][j] / 2 + A3[i][j] / 6
               for j in range(n)] for i in range(n)]
    r = matrix_exp(A)
    assert max(abs(r[i][j] - approx[i][j]) for i in range(n) for j in range(n)) < 1e-4


def test_validation():
    with pytest.raises(ValueError):
        matrix_exp([[1.0, 2.0]])
    with pytest.raises(ValueError):
        matrix_exp([])
