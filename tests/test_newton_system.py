"""Root finding for nonlinear systems: Newton and Broyden."""

import math

import pytest

from quantforge import newton_system, broyden


def _norm(v):
    return sum(c * c for c in v) ** 0.5


def test_circle_line_intersection():
    def f(v):
        x, y = v
        return [x * x + y * y - 1, x - y]
    sol, _ = newton_system(f, [0.5, 0.9])
    assert _norm(f(sol)) < 1e-9
    assert abs(abs(sol[0]) - math.sqrt(0.5)) < 1e-6


def test_transcendental_system():
    def f(v):
        x, y = v
        return [math.exp(x) + y - 1, x + math.exp(y) - 1]
    sol, _ = newton_system(f, [0.5, 0.5])
    assert _norm(sol) < 1e-6            # unique root at (0, 0)


def test_three_variable_system():
    def f(v):
        x, y, z = v
        return [x + y + z - 6, x * x + y * y + z * z - 14,
                x ** 3 + y ** 3 + z ** 3 - 36]
    sol, _ = newton_system(f, [0.5, 1.5, 2.5])   # root is a permutation of (1,2,3)
    assert _norm(f(sol)) < 1e-8


def test_newton_matches_broyden():
    systems = [
        (lambda v: [v[0] ** 2 + v[1] ** 2 - 1, v[0] - v[1]], [0.5, 0.9]),
        (lambda v: [math.exp(v[0]) + v[1] - 1, v[0] + math.exp(v[1]) - 1], [0.5, 0.5]),
    ]
    for f, x0 in systems:
        sn, _ = newton_system(f, x0)
        sb, _ = broyden(f, x0)
        assert _norm([sn[i] - sb[i] for i in range(len(sn))]) < 1e-6


def test_already_at_root():
    def f(v):
        return [v[0] - 1.0, v[1] - 2.0]
    sol, its = newton_system(f, [1.0, 2.0])
    assert its == 0 and _norm(f(sol)) < 1e-12


def test_singular_jacobian_raises():
    with pytest.raises(ValueError):
        newton_system(lambda v: [1.0, 1.0], [0.0, 0.0])
