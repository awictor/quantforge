"""Tests for Fornberg finite-difference weights: known stencils and polynomial exactness."""

import math
import random

import pytest

from quantforge.fd_weights import fd_weights


def _close(a, b, tol=1e-7):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def test_central_stencils():
    w = fd_weights(0, [-1, 0, 1], 2)
    assert all(_close(a, b) for a, b in zip(w[0], [0, 1, 0]))       # value
    assert all(_close(a, b) for a, b in zip(w[1], [-0.5, 0, 0.5]))  # 1st derivative
    assert all(_close(a, b) for a, b in zip(w[2], [1, -2, 1]))      # 2nd derivative


def test_forward_stencil():
    w = fd_weights(0, [0, 1, 2], 1)
    assert all(_close(a, b) for a, b in zip(w[1], [-1.5, 2, -0.5]))


def test_backward_stencil():
    w = fd_weights(0, [-2, -1, 0], 1)
    assert all(_close(a, b) for a, b in zip(w[1], [0.5, -2, 1.5]))


def test_fuzz_polynomial_exactness():
    rng = random.Random(421)
    for _ in range(3000):
        n = rng.randint(1, 6)
        x0 = rng.uniform(-3, 3)
        offs = rng.sample(range(-10, 11), n + 1)
        grid = [o * 0.5 for o in offs]
        m = rng.randint(0, n)
        w = fd_weights(x0, grid, m)
        for k in range(m + 1):
            for deg in range(n + 1):
                if deg < k:
                    exact = 0.0
                else:
                    coef = 1
                    for j in range(k):
                        coef *= (deg - j)
                    exact = coef * (x0 ** (deg - k))
                approx = sum(w[k][i] * (grid[i] ** deg) for i in range(n + 1))
                assert _close(approx, exact, 1e-4)


def test_applied_to_sin():
    grid = [1 + h for h in (-0.02, -0.01, 0, 0.01, 0.02)]
    w = fd_weights(1, grid, 2)
    d1 = sum(w[1][i] * math.sin(grid[i]) for i in range(5))
    d2 = sum(w[2][i] * math.sin(grid[i]) for i in range(5))
    assert _close(d1, math.cos(1), 1e-6)
    assert _close(d2, -math.sin(1), 1e-6)


def test_non_uniform_grid():
    grid = [-2, -0.5, 1, 3]
    w = fd_weights(0, grid, 1)
    # 1st derivative of x^2 at 0 is 0
    assert _close(sum(w[1][i] * g ** 2 for i, g in enumerate(grid)), 0.0, 1e-9)


def test_weights_sum_to_zero_for_derivatives():
    # any derivative weight set sums to 0 (reproduces the derivative of a constant = 0)
    w = fd_weights(0.3, [-1, 0, 1, 2], 3)
    for k in range(1, 4):
        assert _close(sum(w[k]), 0.0, 1e-9)
    assert _close(sum(w[0]), 1.0)  # value weights sum to 1


def test_zeroth_derivative_is_lagrange_interpolation():
    # w[0] are the Lagrange interpolation weights at x0
    grid = [0, 1, 2, 3]
    x0 = 1.5
    w = fd_weights(x0, grid, 0)
    ys = [2 * x ** 3 - x + 1 for x in grid]
    interp = sum(w[0][i] * ys[i] for i in range(4))
    assert _close(interp, 2 * x0 ** 3 - x0 + 1, 1e-9)


def test_too_few_points_raises():
    with pytest.raises(ValueError):
        fd_weights(0, [0, 1], 3)


def test_duplicate_grid_raises():
    with pytest.raises(ValueError):
        fd_weights(0, [1, 1, 2], 1)


def test_negative_deriv_raises():
    with pytest.raises(ValueError):
        fd_weights(0, [0, 1, 2], -1)
