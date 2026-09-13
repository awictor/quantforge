"""Polynomial interpolation: Neville and Newton divided differences."""

import math
import random

import pytest

from quantforge import neville, divided_differences, newton_polynomial


def test_recovers_cubic_exactly():
    p = lambda x: 2 * x ** 3 - 3 * x ** 2 + x - 5
    xs = [0, 1, 2, 3]
    ys = [p(x) for x in xs]
    for xt in (0.5, 1.5, 2.7, -1.0):
        v, _ = neville(xs, ys, xt)
        assert abs(v - p(xt)) < 1e-9


def test_neville_matches_newton():
    p = lambda x: 2 * x ** 3 - 3 * x ** 2 + x - 5
    xs = [0, 1, 2, 3]
    ys = [p(x) for x in xs]
    coef = divided_differences(xs, ys)
    for xt in (0.5, 1.5, 2.7):
        v, _ = neville(xs, ys, xt)
        assert abs(v - newton_polynomial(xs, coef, xt)) < 1e-12


def test_leading_divided_difference_is_leading_coefficient():
    p = lambda x: 2 * x ** 3 - 3 * x ** 2 + x - 5
    xs = [0, 1, 2, 3]
    coef = divided_differences(xs, [p(x) for x in xs])
    assert abs(coef[-1] - 2.0) < 1e-12


def test_richardson_via_neville():
    # Extrapolate the central-difference derivative of sin at h -> 0.
    cd = lambda f, x, h: (f(x + h) - f(x - h)) / (2 * h)
    hs = [0.4, 0.2, 0.1, 0.05]
    vals = [cd(math.sin, 1.0, h) for h in hs]
    v, _ = neville([h * h for h in hs], vals, 0.0)
    assert abs(v - math.cos(1.0)) < 1e-10


def test_error_estimate_honest():
    xs = [0, 0.5, 1, 1.5, 2]
    ys = [math.exp(x) for x in xs]
    v, e = neville(xs, ys, 0.75)
    assert abs(v - math.exp(0.75)) <= 10 * e


def test_random_polynomials_exact():
    rng = random.Random(3)
    for _ in range(200):
        deg = rng.randint(1, 6)
        cf = [rng.uniform(-3, 3) for _ in range(deg + 1)]
        poly = lambda x: sum(cf[k] * x ** k for k in range(deg + 1))
        xs = [i * 0.5 for i in range(deg + 1)]
        ys = [poly(x) for x in xs]
        xt = rng.uniform(-1, 3)
        v, _ = neville(xs, ys, xt)
        assert abs(v - poly(xt)) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        neville([1, 2], [1], 0.5)                 # length mismatch
    with pytest.raises(ValueError):
        neville([1, 1], [2, 3], 0.5)              # duplicate x
    with pytest.raises(ValueError):
        divided_differences([1, 1], [2, 3])        # duplicate x
