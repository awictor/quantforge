import math

from quantforge import (
    gp_predict,
    gp_log_marginal_likelihood,
    rbf_kernel,
    matern32_kernel,
)
from quantforge.linalg import cholesky
from quantforge.gaussian_process import _cho_solve


def close(a, b, tol=1e-4):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _data():
    X = [0.0, 1.0, 2.0, 3.0, 4.0]
    return X, [math.sin(x) for x in X]


def test_interpolates_training_points():
    X, y = _data()
    pred = gp_predict(X, y, X, length_scale=1.0, variance=1.0, noise=1e-8)
    for i in range(len(X)):
        assert close(pred["mean"][i], y[i], 1e-4)
        assert pred["var"][i] < 1e-4


def test_variance_grows_away_from_data():
    X, y = _data()
    p = gp_predict(X, y, [2.0, 10.0], length_scale=1.0, variance=1.0, noise=1e-8)
    assert p["var"][1] > p["var"][0]
    assert close(p["var"][1], 1.0, 0.05)   # -> prior variance far from data


def test_log_marginal_likelihood_matches_brute():
    X, y = _data()
    n = len(X)
    ls, v, noise = 1.0, 1.0, 0.1
    K = [[rbf_kernel(X[i], X[j], ls, v) + (noise if i == j else 0.0) for j in range(n)]
         for i in range(n)]
    M = [row[:] for row in K]
    logdet = 0.0
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[p] = M[p], M[k]
        logdet += math.log(abs(M[k][k]))
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
    L = cholesky(K)
    a = _cho_solve(L, y)
    brute = -0.5 * sum(y[i] * a[i] for i in range(n)) - 0.5 * logdet - 0.5 * n * math.log(2 * math.pi)
    assert close(gp_log_marginal_likelihood(X, y, length_scale=ls, variance=v, noise=noise), brute, 1e-6)


def test_marginal_likelihood_selects_length_scale():
    X, y = _data()
    lmls = {ls: gp_log_marginal_likelihood(X, y, length_scale=ls, variance=1.0, noise=0.01)
            for ls in (0.01, 1.0, 100.0)}
    assert lmls[1.0] > lmls[0.01]
    assert lmls[1.0] > lmls[100.0]


def test_matern_kernel():
    assert close(matern32_kernel(1.0, 1.0, 1.0, 2.0), 2.0)   # r=0 -> variance
    X, y = _data()
    pm = gp_predict(X, y, X, kernel=matern32_kernel, length_scale=1.0, noise=1e-8)
    for i in range(len(X)):
        assert close(pm["mean"][i], y[i], 1e-3)


def test_two_dimensional_inputs():
    X = [[0, 0], [1, 0], [0, 1], [1, 1]]
    y = [0.0, 1.0, 1.0, 2.0]
    p = gp_predict(X, y, [[0.5, 0.5]], length_scale=1.0, noise=1e-6)
    assert 0.0 < p["mean"][0] < 2.0
