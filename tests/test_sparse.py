import random

import pytest

from quantforge import CSRMatrix
from quantforge.gmres import gmres
from quantforge.lu import lu_solve


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_matvec_matches_dense():
    random.seed(0)
    m, n = 8, 6
    A = [[random.gauss(0, 1) if random.random() < 0.4 else 0.0 for _ in range(n)] for _ in range(m)]
    S = CSRMatrix.from_dense(A)
    x = [random.gauss(0, 1) for _ in range(n)]
    dense_mv = [sum(A[i][j] * x[j] for j in range(n)) for i in range(m)]
    for a, c in zip(S.matvec(x), dense_mv):
        assert close(a, c)
    for a, c in zip(S @ x, dense_mv):
        assert close(a, c)


def test_roundtrip_to_dense():
    random.seed(0)
    A = [[random.gauss(0, 1) if random.random() < 0.4 else 0.0 for _ in range(6)] for _ in range(8)]
    D = CSRMatrix.from_dense(A).to_dense()
    for i in range(8):
        for j in range(6):
            assert close(D[i][j], A[i][j])


def test_transpose():
    random.seed(0)
    m, n = 8, 6
    A = [[random.gauss(0, 1) if random.random() < 0.4 else 0.0 for _ in range(n)] for _ in range(m)]
    S = CSRMatrix.from_dense(A)
    T = S.transpose()
    assert T.shape == (n, m)
    Td = T.to_dense()
    for i in range(m):
        for j in range(n):
            assert close(Td[j][i], A[i][j])
    y = [random.gauss(0, 1) for _ in range(m)]
    dense_Tmv = [sum(A[i][j] * y[i] for i in range(m)) for j in range(n)]
    for a, c in zip(T.matvec(y), dense_Tmv):
        assert close(a, c)


def test_triplets_sum_duplicates():
    S = CSRMatrix.from_triplets([(0, 0, 1.0), (0, 0, 2.0), (1, 2, 5.0), (2, 1, -3.0)], (3, 3))
    D = S.to_dense()
    assert close(D[0][0], 3.0)
    assert close(D[1][2], 5.0)
    assert close(D[2][1], -3.0)
    assert S.nnz == 3


def test_gmres_with_csr_operator():
    random.seed(2)
    n = 15
    Ad = [[0.0] * n for _ in range(n)]
    for i in range(n):
        Ad[i][i] = 4.0
        if i > 0:
            Ad[i][i - 1] = Ad[i - 1][i] = -1.0
    Asp = CSRMatrix.from_dense(Ad)
    b = [random.gauss(0, 1) for _ in range(n)]
    res = gmres(Asp.matvec, b, tol=1e-12)
    xref = lu_solve(Ad, b)
    for a, c in zip(res["x"], xref):
        assert close(a, c, 1e-6)


def test_sparse_storage_smaller_than_dense():
    n = 15
    Ad = [[0.0] * n for _ in range(n)]
    for i in range(n):
        Ad[i][i] = 4.0
        if i > 0:
            Ad[i][i - 1] = Ad[i - 1][i] = -1.0
    assert CSRMatrix.from_dense(Ad).nnz < n * n


def test_triplet_out_of_range_raises():
    with pytest.raises(ValueError):
        CSRMatrix.from_triplets([(5, 0, 1.0)], (3, 3))
