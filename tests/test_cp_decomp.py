import random

import pytest

from quantforge import cp_als, cp_reconstruct


def _rand_tensor(I, J, K, R, seed=0):
    random.seed(seed)
    A = [[random.random() for _ in range(R)] for _ in range(I)]
    B = [[random.random() for _ in range(R)] for _ in range(J)]
    C = [[random.random() for _ in range(R)] for _ in range(K)]
    return cp_reconstruct(A, B, C)


def _fro(X):
    return sum(X[i][j][k] ** 2 for i in range(len(X)) for j in range(len(X[0]))
               for k in range(len(X[0][0]))) ** 0.5


def test_recovers_rank_two_tensor():
    X = _rand_tensor(4, 5, 3, 2)
    res = cp_als(X, 2, iterations=500, seed=1)
    assert res["error"] / _fro(X) < 0.02


def test_reconstruction_matches():
    X = _rand_tensor(4, 5, 3, 2)
    res = cp_als(X, 2, iterations=500, seed=1)
    Xh = cp_reconstruct(res["A"], res["B"], res["C"])
    me = max(abs(X[i][j][k] - Xh[i][j][k]) for i in range(4) for j in range(5) for k in range(3))
    assert me / _fro(X) < 0.02


def test_error_decreases():
    X = _rand_tensor(4, 5, 3, 2)
    e10 = cp_als(X, 2, iterations=10, seed=1)["error"]
    e300 = cp_als(X, 2, iterations=300, seed=1)["error"]
    assert e300 <= e10


def test_rank_capacity():
    X = _rand_tensor(4, 5, 3, 3, seed=2)
    e2 = cp_als(X, 2, iterations=300, seed=2)["error"]
    e3 = cp_als(X, 3, iterations=300, seed=2)["error"]
    assert e3 <= e2


def test_reproducible():
    X = _rand_tensor(4, 5, 3, 2)
    a = cp_als(X, 2, iterations=50, seed=7)
    b = cp_als(X, 2, iterations=50, seed=7)
    assert a["A"] == b["A"]


def test_invalid_rank_raises():
    X = _rand_tensor(4, 5, 3, 2)
    with pytest.raises(ValueError):
        cp_als(X, 0)
