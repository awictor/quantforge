import random

import pytest

from quantforge import hankel_singular_values, balanced_truncation
from quantforge.state_space import lti_step_response, dc_gain
from quantforge.lu import lu_solve


def close(a, b, tol=1e-5):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


A = [[0.5, 0.1, 0.0], [0.0, 0.4, 0.05], [0.0, 0.0, 0.3]]
B = [[1.0], [0.5], [0.1]]
C = [[1.0, 0.2, 0.05]]
D = [[0.0]]


def test_hankel_singular_values_descending():
    hsv = hankel_singular_values(A, B, C)
    assert len(hsv) == 3
    assert all(hsv[i] >= hsv[i + 1] >= 0 for i in range(2))


def test_full_order_truncation_preserves_system():
    Ar, Br, Cr, Dr = balanced_truncation(A, B, C, D, 3)
    s_full = [r[0] for r in lti_step_response(A, B, C, D, 30)]
    s_red = [r[0] for r in lti_step_response(Ar, Br, Cr, Dr, 30)]
    assert max(abs(s_full[k] - s_red[k]) for k in range(30)) < 1e-6


def test_dc_gain_preserved_full_order():
    Ar, Br, Cr, Dr = balanced_truncation(A, B, C, D, 3)
    assert close(dc_gain(Ar, Br, Cr, Dr)[0][0], dc_gain(A, B, C, D)[0][0])


def test_reduced_order_keeps_dominant():
    Ar, Br, Cr, Dr = balanced_truncation(A, B, C, D, 2)
    assert len(Ar) == 2
    assert abs(dc_gain(Ar, Br, Cr, Dr)[0][0] - dc_gain(A, B, C, D)[0][0]) < 0.5


def test_hsv_invariant_under_similarity():
    random.seed(0)
    n = 3
    T = [[random.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    cols = [lu_solve(T, [1.0 if i == k else 0.0 for i in range(n)]) for k in range(n)]
    Ti = [[cols[j][i] for j in range(n)] for i in range(n)]

    def mm(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(len(Y))) for j in range(len(Y[0]))]
                for i in range(len(X))]

    At = mm(mm(Ti, A), T)
    Bt = mm(Ti, B)
    Ct = mm(C, T)
    for a, b in zip(hankel_singular_values(A, B, C), hankel_singular_values(At, Bt, Ct)):
        assert close(a, b, 1e-4)


def test_invalid_order_raises():
    with pytest.raises(ValueError):
        balanced_truncation(A, B, C, D, 0)
    with pytest.raises(ValueError):
        balanced_truncation(A, B, C, D, 5)
