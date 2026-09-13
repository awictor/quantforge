"""1-D wave equation (explicit finite difference)."""

import math

import pytest

from quantforge import wave_equation


def _grid(N=201, L=1.0):
    dx = L / (N - 1)
    return N, L, dx, [i * dx for i in range(N)]


def test_standing_wave_matches_dalembert():
    N, L, dx, x = _grid()
    c = 1.0
    u0 = [math.sin(math.pi * xi / L) for xi in x]
    dt = 0.5 * dx / c
    nsteps = int(0.5 / dt)
    u = wave_equation(u0, [0.0] * N, c, dx, dt, nsteps)
    T = nsteps * dt
    exact = [math.cos(math.pi * c * T / L) * math.sin(math.pi * xi / L) for xi in x]
    assert max(abs(u[i] - exact[i]) for i in range(N)) < 1e-3


def test_full_period_returns():
    N, L, dx, x = _grid()
    c = 1.0
    u0 = [math.sin(math.pi * xi / L) for xi in x]
    dt = 0.5 * dx / c
    nsteps = int(2 * L / c / dt)
    u = wave_equation(u0, [0.0] * N, c, dx, dt, nsteps)
    assert max(abs(u[i] - u0[i]) for i in range(N)) < 1e-2


def test_boundaries_fixed():
    N, L, dx, x = _grid()
    u0 = [math.sin(math.pi * xi / L) for xi in x]
    u = wave_equation(u0, [0.0] * N, 1.0, dx, 0.5 * dx, 50)
    assert u[0] == 0.0 and u[-1] == 0.0


def test_pulse_splits_into_halves():
    N, L, dx, x = _grid()
    c = 1.0
    u0 = [math.exp(-((xi - 0.5) ** 2) / 0.002) for xi in x]
    dt = 0.5 * dx / c
    u = wave_equation(u0, [0.0] * N, c, dx, dt, int(0.2 / dt))
    assert abs(max(u) - 0.5) < 0.05        # d'Alembert: two ~0.5 halves


def test_cfl_violation_raises():
    N, L, dx, x = _grid()
    u0 = [0.0] * N
    with pytest.raises(ValueError):
        wave_equation(u0, [0.0] * N, 1.0, dx, 2 * dx, 10)   # C = 2 > 1


def test_validation():
    with pytest.raises(ValueError):
        wave_equation([1.0, 2.0], [0.0, 0.0], 1.0, 0.1, 0.05, 1)       # too few points
    with pytest.raises(ValueError):
        wave_equation([0.0, 1.0, 0.0], [0.0, 0.0], 1.0, 0.1, 0.05, 1)  # v0 length mismatch
