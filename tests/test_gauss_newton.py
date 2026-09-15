import math

import pytest

from quantforge import gauss_newton


def close(a, b, tol=1e-5):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_linear_least_squares_matches_ols():
    xs = [i * 0.3 for i in range(20)]
    ys = [2.5 * x - 1.3 for x in xs]

    def resid(p):
        a, b = p
        return [a * xs[i] + b - ys[i] for i in range(len(xs))]

    r = gauss_newton(resid, [0.0, 0.0])
    assert r["converged"]
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(xs[i] * ys[i] for i in range(n))
    ahat = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    bhat = (sy - ahat * sx) / n
    assert close(r["p"][0], ahat, 1e-6)
    assert close(r["p"][1], bhat, 1e-6)


def test_exponential_decay_fit():
    A_true, k_true = 3.0, 0.8
    ts = [i * 0.25 for i in range(30)]
    ys = [A_true * math.exp(-k_true * t) for t in ts]

    def resid(p):
        A, k = p
        return [A * (-(k * ts[i])).exp() - ys[i] for i in range(len(ts))]

    r = gauss_newton(resid, [1.0, 0.1])
    assert r["converged"]
    assert close(r["p"][0], A_true, 1e-4)
    assert close(r["p"][1], k_true, 1e-4)
    assert r["cost"] < 1e-12


def test_rosenbrock_as_residuals():
    def resid(p):
        x, y = p
        return [1 - x, 10 * (y - x * x)]

    r = gauss_newton(resid, [-1.2, 1.0])
    assert close(r["p"][0], 1.0, 1e-6)
    assert close(r["p"][1], 1.0, 1e-6)
    assert r["cost"] < 1e-16


def test_sinusoid_fit():
    At, wt, pt = 2.0, 1.5, 0.5
    tt = [i * 0.2 for i in range(40)]
    yy = [At * math.sin(wt * t + pt) for t in tt]

    def resid(p):
        A, w, ph = p
        return [A * (w * tt[i] + ph).sin() - yy[i] for i in range(len(tt))]

    r = gauss_newton(resid, [1.5, 1.4, 0.4])
    assert r["cost"] < 1e-10


def test_empty_parameters_raises():
    with pytest.raises(ValueError):
        gauss_newton(lambda p: [p[0]], [])
