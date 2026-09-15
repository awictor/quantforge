import random

import pytest

from quantforge import nmf
from quantforge.nmf import _matmul


def _fro(A, B):
    return sum((A[i][j] - B[i][j]) ** 2 for i in range(len(A)) for j in range(len(A[0]))) ** 0.5


def _low_rank_nonneg(m, n, r, seed=0):
    random.seed(seed)
    W = [[random.random() for _ in range(r)] for _ in range(m)]
    H = [[random.random() for _ in range(n)] for _ in range(r)]
    return _matmul(W, H)


def test_recovers_low_rank_nonneg_product():
    V = _low_rank_nonneg(15, 10, 3)
    res = nmf(V, 3, iterations=2000, seed=1)
    vnorm = sum(V[i][j] ** 2 for i in range(15) for j in range(10)) ** 0.5
    assert res["error"] / vnorm < 0.02


def test_factors_are_nonnegative():
    V = _low_rank_nonneg(15, 10, 3)
    res = nmf(V, 3, iterations=500, seed=1)
    assert all(w >= 0 for row in res["W"] for w in row)
    assert all(h >= 0 for row in res["H"] for h in row)


def test_reconstruction_matches():
    V = _low_rank_nonneg(15, 10, 3)
    res = nmf(V, 3, iterations=2000, seed=1)
    vnorm = sum(V[i][j] ** 2 for i in range(15) for j in range(10)) ** 0.5
    assert _fro(V, _matmul(res["W"], res["H"])) / vnorm < 0.02


def test_error_decreases_with_iterations():
    V = _low_rank_nonneg(15, 10, 3)
    e50 = nmf(V, 3, iterations=50, seed=1)["error"]
    e500 = nmf(V, 3, iterations=500, seed=1)["error"]
    assert e500 <= e50


def test_higher_k_lowers_error():
    V = _low_rank_nonneg(15, 10, 3)
    e2 = nmf(V, 2, iterations=500, seed=1)["error"]
    e4 = nmf(V, 4, iterations=500, seed=1)["error"]
    assert e4 <= e2


def test_reproducible():
    V = _low_rank_nonneg(15, 10, 3)
    a = nmf(V, 3, iterations=100, seed=7)
    b = nmf(V, 3, iterations=100, seed=7)
    assert a["W"] == b["W"] and a["H"] == b["H"]


def test_block_structure_recovered():
    V = [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]]
    res = nmf(V, 2, iterations=2000, seed=3)
    WH = _matmul(res["W"], res["H"])
    assert max(abs(V[i][j] - WH[i][j]) for i in range(4) for j in range(4)) < 0.05


def test_negative_input_raises():
    with pytest.raises(ValueError):
        nmf([[-1, 2], [3, 4]], 2)
