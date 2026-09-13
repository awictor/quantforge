"""1-D heat/diffusion equation (Crank-Nicolson)."""

import math

import pytest

from quantforge import heat_equation_cn


def test_sine_mode_decays_analytically():
    L, N = 1.0, 101
    dx = L / (N - 1)
    alpha = 0.5
    x = [i * dx for i in range(N)]
    u0 = [math.sin(math.pi * xi / L) for xi in x]
    dt, nsteps = 0.0005, 200
    T = nsteps * dt
    u = heat_equation_cn(u0, alpha, dx, dt, nsteps, left=0.0, right=0.0)
    decay = math.exp(-alpha * (math.pi / L) ** 2 * T)
    assert max(abs(u[i] - decay * math.sin(math.pi * x[i] / L)) for i in range(N)) < 1e-3


def test_boundaries_held():
    u0 = [0.0] * 50
    u = heat_equation_cn(u0, 0.5, 0.02, 0.01, 10, left=1.0, right=2.0)
    assert u[0] == 1.0 and u[-1] == 2.0


def test_steady_state_is_linear():
    N = 101
    L = 1.0
    dx = L / (N - 1)
    x = [i * dx for i in range(N)]
    u = heat_equation_cn([0.0] * N, 0.5, dx, 0.01, 5000, left=0.0, right=10.0)
    lin = [10.0 * xi / L for xi in x]
    assert max(abs(u[i] - lin[i]) for i in range(N)) < 1e-3


def test_bump_spreads():
    N = 101
    dx = 1.0 / (N - 1)
    x = [i * dx for i in range(N)]
    u0 = [math.exp(-((xi - 0.5) ** 2) / 0.01) for xi in x]
    u = heat_equation_cn(u0, 0.5, dx, 0.001, 100, left=0.0, right=0.0)
    assert max(u) < max(u0)                # diffusion lowers the peak


def test_unconditional_stability_large_dt():
    N = 101
    dx = 1.0 / (N - 1)
    x = [i * dx for i in range(N)]
    u0 = [math.exp(-((xi - 0.5) ** 2) / 0.01) for xi in x]
    u = heat_equation_cn(u0, 0.5, dx, 1.0, 10, left=0.0, right=0.0)
    assert all(abs(v) < 10 for v in u)


def test_validation():
    with pytest.raises(ValueError):
        heat_equation_cn([1.0, 2.0], 0.5, 0.1, 0.1, 1)
    with pytest.raises(ValueError):
        heat_equation_cn([0.0, 1.0, 0.0], 0.5, 0.0, 0.1, 1)
