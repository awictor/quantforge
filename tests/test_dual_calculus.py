"""Multivariate autodiff: exact gradients and autodiff Newton."""

import math

import pytest

from quantforge import dual_gradient, dual_newton, gradient
from quantforge.dual import exp, log, sqrt
from quantforge.rootfind import newton


def test_gradient_exact():
    f = lambda v: v[0] * v[0] + 3 * v[1] * v[1] + v[0] * v[1]
    g = dual_gradient(f, [1.0, 2.0])
    assert abs(g[0] - 4.0) < 1e-12 and abs(g[1] - 13.0) < 1e-12


def test_gradient_matches_numerical():
    f = lambda v: v[0] * v[0] + 3 * v[1] * v[1] + v[0] * v[1]
    fp = lambda v: v[0] ** 2 + 3 * v[1] ** 2 + v[0] * v[1]
    g = dual_gradient(f, [1.0, 2.0])
    gn = gradient(fp, [1.0, 2.0])
    assert max(abs(g[i] - gn[i]) for i in range(2)) < 1e-6


def test_gradient_transcendental():
    f = lambda v: exp(v[0]) * log(v[1]) + sqrt(v[2])
    g = dual_gradient(f, [0.5, 3.0, 4.0])
    ana = [math.exp(0.5) * math.log(3.0), math.exp(0.5) / 3.0, 0.5 / math.sqrt(4.0)]
    assert max(abs(g[i] - ana[i]) for i in range(3)) < 1e-12


def test_newton_sqrt2():
    r = dual_newton(lambda x: x * x - 2, 1.0)
    assert r["converged"]
    assert abs(r["root"] - math.sqrt(2)) < 1e-12


def test_newton_exp():
    r = dual_newton(lambda x: exp(x) - 2, 0.0)
    assert abs(r["root"] - math.log(2)) < 1e-12


def test_newton_matches_library():
    d = dual_newton(lambda x: x * x - 2, 1.0)["root"]
    lib = newton(lambda x: x * x - 2, lambda x: 2 * x, 1.0)
    assert abs(d - lib) < 1e-12


def test_newton_zero_derivative_raises():
    with pytest.raises(ValueError):
        dual_newton(lambda x: x * 0 + 1, 1.0)     # constant -> zero derivative
