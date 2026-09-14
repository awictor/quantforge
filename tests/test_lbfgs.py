import math
import random

import pytest

from quantforge import lbfgs, bfgs, reverse_gradient


def close(a, b, tol=1e-4):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_quadratic_bowl():
    target = [3.0, -2.0, 5.0]

    def quad(x):
        return sum((xi - ti) ** 2 for xi, ti in zip(x, target)) + 1.0

    r = lbfgs(quad, [0.0, 0.0, 0.0])
    assert r["converged"]
    for xi, ti in zip(r["x"], target):
        assert close(xi, ti)
    assert close(r["fun"], 1.0)


def test_rosenbrock():
    def rosen(x):
        return (1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2

    r = lbfgs(rosen, [-1.2, 1.0], max_iter=2000)
    assert close(r["x"][0], 1.0, 1e-3)
    assert close(r["x"][1], 1.0, 1e-3)
    assert r["fun"] < 1e-8


def test_agrees_with_bfgs():
    def f(x):
        return (x[0] - 1) ** 2 + (x[1] + 2) ** 2 + 0.5 * x[0] * x[1]

    rl = lbfgs(f, [5.0, 5.0])
    rb = bfgs(f, [5.0, 5.0])
    for a, b in zip(rl["x"], rb["x"]):
        assert close(a, b, 1e-3)


def test_analytic_gradient_matches_finite_difference():
    def fr(v):
        return (v[0] - 1) ** 2 + (v[1] + 2) ** 2 + 0.5 * v[0] * v[1]

    r_an = lbfgs(fr, [5.0, 5.0], grad=lambda x: reverse_gradient(fr, x))
    r_fd = lbfgs(fr, [5.0, 5.0])
    assert r_an["converged"]
    for a, b in zip(r_an["x"], r_fd["x"]):
        assert close(a, b, 1e-4)


def test_high_dimensional_quadratic():
    random.seed(0)
    t = [random.uniform(-5, 5) for _ in range(100)]

    def q100(x):
        return sum((xi - ti) ** 2 for xi, ti in zip(x, t))

    r = lbfgs(q100, [0.0] * 100, m=8)
    assert r["converged"]
    assert max(abs(a - b) for a, b in zip(r["x"], t)) < 1e-3


def test_strong_wolfe_on_mixed_objective():
    def g(x):
        return math.exp(0.3 * x[0]) + x[0] ** 2 + (x[1] - 3) ** 2

    r = lbfgs(g, [2.0, 0.0])
    assert r["converged"]
    h = 1e-6
    gx = (g([r["x"][0] + h, r["x"][1]]) - g([r["x"][0] - h, r["x"][1]])) / (2 * h)
    assert abs(gx) < 1e-4


def test_empty_dimension_raises():
    with pytest.raises(ValueError):
        lbfgs(lambda x: 0.0, [])
