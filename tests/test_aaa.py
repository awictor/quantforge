import math

import pytest

from quantforge import aaa


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_recovers_true_rational():
    f = lambda x: (x + 2) / (x * x + 1)
    xs = [-3 + 6 * i / 100 for i in range(101)]
    ys = [f(x) for x in xs]
    r = aaa(xs, ys, tol=1e-12)
    assert max(abs(r(x) - f(x)) for x in (-2.5, -1.3, 0.2, 1.7, 2.9)) < 1e-8
    assert len(r.support_x) <= 5


def test_interpolates_support_points():
    f = lambda x: (x + 2) / (x * x + 1)
    xs = [-3 + 6 * i / 100 for i in range(101)]
    r = aaa(xs, [f(x) for x in xs], tol=1e-12)
    for i, sx in enumerate(r.support_x):
        assert close(r(sx), r.support_y[i])


def test_approximates_exp():
    xs = [-1 + 2 * i / 200 for i in range(201)]
    r = aaa(xs, [math.exp(x) for x in xs], tol=1e-10)
    assert max(abs(r(x) - math.exp(x)) for x in (-0.9, -0.3, 0.5, 0.95)) < 1e-8


def test_near_pole_and_capture():
    f = lambda x: 1 / (1.05 - x)
    xs = [-1 + 2 * i / 300 for i in range(301)]
    r = aaa(xs, [f(x) for x in xs], tol=1e-10)
    assert max(abs(r(x) - f(x)) for x in (-0.8, 0.0, 0.7, 0.99)) < 1e-6
    assert abs(r(1.04)) > 10          # captures the pole just outside the sample range


def test_approximates_log():
    xs = [0.1 + 2.9 * i / 200 for i in range(201)]
    r = aaa(xs, [math.log(x) for x in xs], tol=1e-10)
    assert max(abs(r(x) - math.log(x)) for x in (0.5, 1.5, 2.5)) < 1e-7


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        aaa([1.0, 2.0], [1.0])
