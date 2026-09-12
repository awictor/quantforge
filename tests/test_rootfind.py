"""General-purpose scalar root finders."""

import math

import pytest

from quantforge import bisection, brent, newton


def test_bisection_sqrt2():
    assert bisection(lambda x: x * x - 2, 0, 2) == pytest.approx(math.sqrt(2), abs=1e-9)


def test_brent_sqrt2():
    assert brent(lambda x: x * x - 2, 0, 2) == pytest.approx(math.sqrt(2), abs=1e-10)


def test_newton_sqrt2():
    assert newton(lambda x: x * x - 2, lambda x: 2 * x, 1.0) == pytest.approx(
        math.sqrt(2), abs=1e-10)


def test_newton_with_bracket():
    assert newton(lambda x: x * x - 2, lambda x: 2 * x, 1.0, lo=0, hi=2) == \
        pytest.approx(math.sqrt(2), abs=1e-10)


def test_transcendental_cos_fixed_point():
    r = 0.7390851332151607
    g = lambda x: math.cos(x) - x
    assert brent(g, 0, 1) == pytest.approx(r, abs=1e-10)
    assert newton(g, lambda x: -math.sin(x) - 1, 0.5) == pytest.approx(r, abs=1e-10)


def test_brent_cubic():
    assert brent(lambda x: x ** 3 - x - 2, 1, 2) == pytest.approx(1.5213797, abs=1e-6)


def test_bisection_endpoint_root():
    assert bisection(lambda x: x, -1, 0) == 0.0


def test_validation():
    with pytest.raises(ValueError):
        brent(lambda x: x * x - 2, 3, 4)      # no sign change
    with pytest.raises(ValueError):
        newton(lambda x: x * x + 1, lambda x: 2 * x, 0.0)   # zero derivative, no bracket
