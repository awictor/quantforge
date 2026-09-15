import math
import random

from quantforge import randomized_svd
from quantforge.randomized_svd import _matmul, _transpose
from quantforge.svd import svd


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _low_rank(m, n, r, seed=0):
    random.seed(seed)
    L = [[random.gauss(0, 1) for _ in range(r)] for _ in range(m)]
    R = [[random.gauss(0, 1) for _ in range(n)] for _ in range(r)]
    return _matmul(L, R)


def test_exact_reconstruction_at_true_rank():
    A = _low_rank(20, 15, 3)
    U, s, V = randomized_svd(A, 3, n_power=3, seed=1)
    D = [[s[i] if i == j else 0.0 for j in range(3)] for i in range(3)]
    Ahat = _matmul(_matmul(U, D), _transpose(V))
    maxerr = max(abs(A[i][j] - Ahat[i][j]) for i in range(20) for j in range(15))
    assert maxerr < 1e-6


def test_singular_values_match_full_svd():
    A = _low_rank(20, 15, 3)
    _, s, _ = randomized_svd(A, 3, n_power=3, seed=1)
    _, sf, _ = svd(A)
    for a, b in zip(sorted(s, reverse=True), sorted(sf, reverse=True)[:3]):
        assert close(a, b)


def test_orthonormal_factors():
    A = _low_rank(20, 15, 3)
    U, s, V = randomized_svd(A, 3, n_power=3, seed=1)
    for j in range(3):
        assert close(sum(U[i][j] ** 2 for i in range(20)), 1.0)
        assert close(sum(V[i][j] ** 2 for i in range(15)), 1.0)


def test_low_rank_approx_of_full_rank_matrix():
    random.seed(2)
    B = [[random.gauss(0, 1) + 5 * math.exp(-0.5 * abs(i - j)) for j in range(12)]
         for i in range(18)]
    _, s, _ = randomized_svd(B, 5, n_power=4, seed=3)
    _, sf, _ = svd(B)
    for a, b in zip(sorted(s, reverse=True), sorted(sf, reverse=True)[:5]):
        assert close(a, b, 1e-2)


def test_reproducible():
    A = _low_rank(20, 15, 3)
    _, sa, _ = randomized_svd(A, 3, seed=7)
    _, sb, _ = randomized_svd(A, 3, seed=7)
    assert sa == sb
