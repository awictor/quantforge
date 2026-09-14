"""One-dimensional SDE integrators: Euler-Maruyama and Milstein."""

import math

import pytest

from quantforge import euler_maruyama, milstein, gbm_paths
from quantforge.sde import _normals


def test_zero_diffusion_is_deterministic():
    drift = lambda x, t: 0.5 * x
    diff = lambda x, t: 0.0
    p = euler_maruyama(drift, diff, 1.0, 1.0, 10000)
    assert abs(p[-1] - math.exp(0.5)) < 1e-3


def test_gbm_terminal_mean():
    paths = gbm_paths(0.1, 0.2, 100.0, 1.0, 50, 3000, seed=42)
    terminal = [p[-1] for p in paths]
    mean = sum(terminal) / len(terminal)
    assert abs(mean - 100 * math.exp(0.1)) < 3


def test_milstein_beats_euler_strong_error():
    mu, sigma, x0, T, nsteps = 0.1, 0.3, 1.0, 1.0, 50
    dt = T / nsteps
    sdt = math.sqrt(dt)
    z = _normals(7, nsteps)
    W = 0.0
    xe = xm = x0
    tc = 0.0
    for k in range(nsteps):
        dw = sdt * z[k]
        xe = xe + mu * xe * dt + sigma * xe * dw
        xm = xm + mu * xm * dt + sigma * xm * dw + 0.5 * sigma * xm * sigma * (dw * dw - dt)
        W += dw
        tc += dt
    exact = x0 * math.exp((mu - 0.5 * sigma * sigma) * tc + sigma * W)
    assert abs(xm - exact) < abs(xe - exact)


def test_reproducible():
    a = gbm_paths(0.1, 0.2, 100, 1, 20, 5, seed=1)
    b = gbm_paths(0.1, 0.2, 100, 1, 20, 5, seed=1)
    assert a == b


def test_validation():
    with pytest.raises(ValueError):
        euler_maruyama(lambda x, t: x, lambda x, t: 0.0, 1, 1, 0)
    with pytest.raises(ValueError):
        gbm_paths(0.1, 0.2, 100, 1, 10, 2, scheme="bad")
