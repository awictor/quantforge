"""BFGS quasi-Newton minimization."""

import math
import random

import pytest

from quantforge import bfgs
from quantforge.portopt import _invert


def test_quadratic_bowl_fast():
    r = bfgs(lambda v: (v[0] - 3) ** 2 + (v[1] + 1) ** 2, [0.0, 0.0])
    assert r["converged"]
    assert abs(r["x"][0] - 3) < 1e-5 and abs(r["x"][1] + 1) < 1e-5
    assert r["n_iter"] <= 5              # superlinear on a quadratic


def test_rosenbrock():
    r = bfgs(lambda v: (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2,
             [-1.2, 1.0], max_iter=1000)
    assert abs(r["x"][0] - 1) < 1e-4 and abs(r["x"][1] - 1) < 1e-4


def test_5d_quadratic_matches_linear_solve():
    rng = random.Random(3)
    A = [[rng.gauss(0, 1) for _ in range(5)] for _ in range(5)]
    M = [[sum(A[k][i] * A[k][j] for k in range(5)) + (1 if i == j else 0)
          for j in range(5)] for i in range(5)]
    b = [rng.gauss(0, 1) for _ in range(5)]
    f = (lambda v: sum(0.5 * v[i] * M[i][j] * v[j] for i in range(5) for j in range(5))
         - sum(b[i] * v[i] for i in range(5)))
    r = bfgs(f, [0.0] * 5, max_iter=500)
    inv = _invert(M)
    xstar = [sum(inv[i][j] * b[j] for j in range(5)) for i in range(5)]
    assert max(abs(r["x"][i] - xstar[i]) for i in range(5)) < 1e-6


def test_convex_smooth():
    lse = lambda v: math.log(sum(math.exp(x) for x in v))
    r = bfgs(lambda v: lse(v) + 0.5 * sum(x * x for x in v), [1.0, -1.0, 0.5])
    assert r["converged"]
    assert r["grad_norm"] < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        bfgs(lambda v: v[0], [])
