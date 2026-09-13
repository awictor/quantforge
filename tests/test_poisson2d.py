"""2-D Poisson equation solver (SOR)."""

import math

import pytest

from quantforge import poisson2d


def test_laplace_linear_boundary_exact():
    N = 21
    dx = dy = 1.0 / (N - 1)
    x = [i * dx for i in range(N)]
    bnd = [[x[i] for i in range(N)] for _ in range(N)]
    for j in range(1, N - 1):
        for i in range(1, N - 1):
            bnd[j][i] = 0.0
    u, _ = poisson2d([[0.0] * N for _ in range(N)], bnd, dx, dy)
    assert max(abs(u[j][i] - x[i]) for j in range(N) for i in range(N)) < 1e-6


def test_poisson_sine_source():
    N = 21
    dx = dy = 1.0 / (N - 1)
    x = [i * dx for i in range(N)]
    ue = lambda xi, yj: math.sin(math.pi * xi) * math.sin(math.pi * yj)
    f = [[-2 * math.pi ** 2 * ue(x[i], x[j]) for i in range(N)] for j in range(N)]
    u, _ = poisson2d(f, [[0.0] * N for _ in range(N)], dx, dy, tol=1e-10, max_iter=20000)
    err = max(abs(u[j][i] - ue(x[i], x[j])) for j in range(N) for i in range(N))
    assert err < 0.01                       # O(h^2) discretization


def test_maximum_principle():
    N = 21
    dx = dy = 1.0 / (N - 1)
    bnd = [[0.0] * N for _ in range(N)]
    for i in range(N):
        bnd[0][i] = 1.0
        bnd[N - 1][i] = 2.0
    for j in range(N):
        bnd[j][0] = 0.5
        bnd[j][N - 1] = 1.5
    u, _ = poisson2d([[0.0] * N for _ in range(N)], bnd, dx, dy)
    interior = [u[j][i] for j in range(1, N - 1) for i in range(1, N - 1)]
    assert min(interior) >= 0.5 - 1e-6 and max(interior) <= 2.0 + 1e-6


def test_sor_faster_than_gauss_seidel():
    N = 21
    dx = dy = 1.0 / (N - 1)
    x = [i * dx for i in range(N)]
    ue = lambda xi, yj: math.sin(math.pi * xi) * math.sin(math.pi * yj)
    f = [[-2 * math.pi ** 2 * ue(x[i], x[j]) for i in range(N)] for j in range(N)]
    _, it_gs = poisson2d(f, [[0.0] * N for _ in range(N)], dx, dy, omega=1.0, max_iter=20000)
    _, it_sor = poisson2d(f, [[0.0] * N for _ in range(N)], dx, dy, omega=1.8, max_iter=20000)
    assert it_sor < it_gs


def test_validation():
    with pytest.raises(ValueError):
        poisson2d([[0, 0], [0, 0]], [[0, 0], [0, 0]], 0.1, 0.1)   # too small
    with pytest.raises(ValueError):
        poisson2d([[0] * 3 for _ in range(3)], [[0] * 3 for _ in range(3)], 0.1, 0.1, omega=2.5)
