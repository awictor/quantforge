import random

import pytest

from quantforge import cur, column_leverage_scores
from quantforge.nmf import _matmul


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _rand(m, n, seed=0):
    random.seed(seed)
    return [[random.gauss(0, 1) for _ in range(n)] for _ in range(m)]


def _low_rank(m, n, r, seed=0):
    random.seed(seed)
    L = [[random.gauss(0, 1) for _ in range(r)] for _ in range(m)]
    R = [[random.gauss(0, 1) for _ in range(n)] for _ in range(r)]
    return _matmul(L, R)


def test_leverage_scores_sum_to_one():
    A = _rand(12, 8)
    lev = column_leverage_scores(A, 3)
    assert close(sum(lev), 1.0, 1e-9)
    assert all(l >= 0 for l in lev)


def test_c_and_r_are_actual_columns_and_rows():
    A = _rand(12, 8)
    res = cur(A, 4, 5, seed=1)
    for jj, j in enumerate(res["col_indices"]):
        for i in range(12):
            assert A[i][j] == res["C"][i][jj]
    for ii, i in enumerate(res["row_indices"]):
        assert A[i] == res["R"][ii]


def test_low_rank_reconstruction():
    V = _low_rank(15, 10, 3)
    vnorm = sum(V[i][j] ** 2 for i in range(15) for j in range(10)) ** 0.5
    res = cur(V, 5, 5, k=3, seed=2)
    assert res["error"] / vnorm < 0.02


def test_more_columns_lower_error():
    V = _low_rank(15, 10, 3)
    e_small = cur(V, 3, 3, k=3, seed=5)["error"]
    e_big = cur(V, 8, 8, k=3, seed=5)["error"]
    assert e_big <= e_small + 1e-9


def test_reproducible():
    A = _rand(12, 8)
    a = cur(A, 4, 5, seed=9)
    b = cur(A, 4, 5, seed=9)
    assert a["col_indices"] == b["col_indices"]
    assert a["row_indices"] == b["row_indices"]


def test_errors():
    A = _rand(12, 8)
    with pytest.raises(ValueError):
        cur(A, 0, 5)
    with pytest.raises(ValueError):
        cur(A, 4, 50)
