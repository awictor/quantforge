import random

from quantforge import tucker_hosvd, tucker_reconstruct
from quantforge.qr import qr_decomposition


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _rand(I, J, K, seed=0):
    random.seed(seed)
    return [[[random.gauss(0, 1) for _ in range(K)] for _ in range(J)] for _ in range(I)]


def _orth(m, r, seed):
    random.seed(seed)
    A = [[random.gauss(0, 1) for _ in range(r)] for _ in range(m)]
    Q, _ = qr_decomposition(A)
    return [[Q[i][j] for j in range(r)] for i in range(m)]


def test_full_rank_exact_reconstruction():
    X = _rand(4, 5, 3)
    res = tucker_hosvd(X)
    Xh = tucker_reconstruct(res["U"], res["V"], res["W"], res["core"])
    me = max(abs(X[i][j][k] - Xh[i][j][k]) for i in range(4) for j in range(5) for k in range(3))
    assert me < 1e-9


def test_factors_orthonormal():
    res = tucker_hosvd(_rand(4, 5, 3))
    for M in (res["U"], res["V"], res["W"]):
        r = len(M[0])
        for a in range(r):
            assert close(sum(M[i][a] ** 2 for i in range(len(M))), 1.0)
            for b in range(a):
                assert abs(sum(M[i][a] * M[i][b] for i in range(len(M)))) < 1e-9


def test_core_dimensions():
    res = tucker_hosvd(_rand(4, 5, 3))
    assert len(res["core"]) == 4
    assert len(res["core"][0]) == 5
    assert len(res["core"][0][0]) == 3


def test_low_multilinear_rank_recovered_by_truncation():
    U0, V0, W0 = _orth(6, 2, 1), _orth(6, 2, 2), _orth(6, 2, 3)
    core = [[[random.gauss(0, 1) for _ in range(2)] for _ in range(2)] for _ in range(2)]
    X = tucker_reconstruct(U0, V0, W0, core)
    res = tucker_hosvd(X, ranks=(2, 2, 2))
    Xh = tucker_reconstruct(res["U"], res["V"], res["W"], res["core"])
    me = max(abs(X[i][j][k] - Xh[i][j][k]) for i in range(6) for j in range(6) for k in range(6))
    assert me < 1e-8


def test_truncation_captures_dominant_energy():
    U0, V0, W0 = _orth(6, 2, 1), _orth(6, 2, 2), _orth(6, 2, 3)
    core = [[[5.0 if a == b == c else 0.0 for c in range(2)] for b in range(2)] for a in range(2)]
    base = tucker_reconstruct(U0, V0, W0, core)
    random.seed(9)
    X = [[[base[i][j][k] + random.gauss(0, 0.05) for k in range(6)] for j in range(6)]
         for i in range(6)]
    res = tucker_hosvd(X, ranks=(2, 2, 2))
    Xh = tucker_reconstruct(res["U"], res["V"], res["W"], res["core"])
    xn = sum(X[i][j][k] ** 2 for i in range(6) for j in range(6) for k in range(6)) ** 0.5
    err = sum((X[i][j][k] - Xh[i][j][k]) ** 2 for i in range(6) for j in range(6)
              for k in range(6)) ** 0.5
    assert err / xn < 0.2
