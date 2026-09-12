"""Numerical differentiation: gradient, Hessian, Jacobian."""

import math

import pytest

from quantforge import gradient, hessian, jacobian


def test_gradient_of_quadratic():
    # f = x0^2 + 3 x1^2 + x0 x1 -> grad = [2x0 + x1, 6x1 + x0]
    f = lambda v: v[0] ** 2 + 3 * v[1] ** 2 + v[0] * v[1]
    g = gradient(f, [1.0, 2.0])
    assert abs(g[0] - 4) < 1e-4
    assert abs(g[1] - 13) < 1e-4


def test_hessian_of_quadratic():
    f = lambda v: v[0] ** 2 + 3 * v[1] ** 2 + v[0] * v[1]
    H = hessian(f, [1.0, 2.0])
    assert abs(H[0][0] - 2) < 1e-2
    assert abs(H[1][1] - 6) < 1e-2
    assert abs(H[0][1] - 1) < 1e-2
    assert abs(H[0][1] - H[1][0]) < 1e-12


def test_gradient_zero_at_minimum():
    g = gradient(lambda v: (v[0] - 1) ** 2 + (v[1] - 2) ** 2, [1.0, 2.0])
    assert all(abs(x) < 1e-6 for x in g)


def test_jacobian_of_vector_map():
    # f = [x0 + x1, x0 * x1] at (2, 3) -> [[1, 1], [3, 2]]
    J = jacobian(lambda v: [v[0] + v[1], v[0] * v[1]], [2.0, 3.0])
    assert abs(J[0][0] - 1) < 1e-4 and abs(J[0][1] - 1) < 1e-4
    assert abs(J[1][0] - 3) < 1e-4 and abs(J[1][1] - 2) < 1e-4


def test_trig_gradient():
    g = gradient(lambda v: math.sin(v[0]) * math.cos(v[1]), [0.0, 0.0])
    assert abs(g[0] - 1) < 1e-4
    assert abs(g[1]) < 1e-6
