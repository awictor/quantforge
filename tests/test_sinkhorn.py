import math
import random

import pytest

from quantforge import sinkhorn, cost_matrix
from quantforge.wasserstein import wasserstein_distance  # noqa: F401  (documents the relation)


def close(a, b, tol=1e-3):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_marginals_and_nonnegativity():
    a = [0.3, 0.5, 0.2]
    b = [0.4, 0.6]
    C = [[1, 2], [3, 1], [2, 2]]
    P = sinkhorn(a, b, C, eps=0.05)["plan"]
    for i in range(3):
        assert close(sum(P[i]), a[i], 1e-5)
    for j in range(2):
        assert close(sum(P[i][j] for i in range(3)), b[j], 1e-5)
    assert all(P[i][j] >= 0 for i in range(3) for j in range(2))


def test_known_transport_small_eps():
    xs = [0.0, 1.0]
    ys = [2.0]
    C = cost_matrix(xs, ys, p=2)
    r = sinkhorn([0.5, 0.5], [1.0], C, eps=0.002, max_iter=20000)
    exact = 0.5 * (2 - 0) ** 2 + 0.5 * (2 - 1) ** 2
    assert close(r["cost"], exact, 2e-2)


def test_matches_1d_wasserstein_for_equal_weight_samples():
    random.seed(0)
    n = 8
    xs = sorted(random.uniform(0, 5) for _ in range(n))
    ys = sorted(random.uniform(0, 5) for _ in range(n))
    C = cost_matrix(xs, ys, p=1)
    r = sinkhorn([1.0 / n] * n, [1.0 / n] * n, C, eps=0.005, max_iter=30000)
    w1 = sum(abs(xs[i] - ys[i]) for i in range(n)) / n
    assert close(r["cost"], w1, 3e-2)


def test_identical_distributions_zero_cost():
    xs = [0.0, 1.0, 2.0]
    a = [1 / 3] * 3
    C = cost_matrix(xs, xs, p=2)
    r = sinkhorn(a, a, C, eps=0.005, max_iter=20000)
    assert r["cost"] < 0.02


def test_symmetry():
    a = [0.2, 0.3, 0.5]
    b = [0.6, 0.4]
    C = [[1.0, 2.5], [2.0, 1.0], [0.5, 3.0]]
    r1 = sinkhorn(a, b, C, eps=0.01)
    CT = [[C[i][j] for i in range(3)] for j in range(2)]
    r2 = sinkhorn(b, a, CT, eps=0.01)
    assert close(r1["cost"], r2["cost"], 1e-4)


def test_vector_cost_matrix():
    xs = [[0.0, 0.0], [1.0, 1.0]]
    ys = [[0.0, 0.0], [1.0, 1.0]]
    C = cost_matrix(xs, ys, p=2)
    assert close(C[0][0], 0.0)
    assert close(C[0][1], 2.0)   # squared distance of (1,1) from origin
    r = sinkhorn([0.5, 0.5], [0.5, 0.5], C, eps=0.01, max_iter=20000)
    assert r["cost"] < 0.05


def test_unequal_mass_raises():
    with pytest.raises(ValueError):
        sinkhorn([0.5, 0.5], [1.0, 1.0], [[1, 1], [1, 1]])
