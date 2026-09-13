"""Simulated annealing global optimization."""

import math

import pytest

from quantforge import simulated_annealing, nelder_mead


def test_sphere():
    r = simulated_annealing(lambda v: sum(x * x for x in v), [3.0, -2.0, 1.0],
                            bounds=[(-5, 5)] * 3, max_iter=20000)
    assert r["fun"] < 1e-3


def test_rosenbrock():
    rosen = lambda v: (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2
    r = simulated_annealing(rosen, [-3.0, 3.0], bounds=[(-5, 5), (-5, 5)],
                            T0=2.0, max_iter=50000)
    assert r["fun"] < 1e-2


def test_escapes_local_minimum():
    rastrigin = lambda v: 10 * len(v) + sum(x * x - 10 * math.cos(2 * math.pi * x) for x in v)
    sa = simulated_annealing(rastrigin, [3.0, 3.0], bounds=[(-5.12, 5.12)] * 2,
                             T0=5.0, max_iter=40000)
    _, nm_f = nelder_mead(rastrigin, [3.0, 3.0], step=0.5)
    assert sa["fun"] < nm_f              # SA climbs out where NM stalls
    assert sa["fun"] < 3.0               # near the global basin


def test_deterministic_for_seed():
    rastrigin = lambda v: 10 * len(v) + sum(x * x - 10 * math.cos(2 * math.pi * x) for x in v)
    r1 = simulated_annealing(rastrigin, [3.0, 3.0], bounds=[(-5.12, 5.12)] * 2,
                             seed=99, max_iter=5000)
    r2 = simulated_annealing(rastrigin, [3.0, 3.0], bounds=[(-5.12, 5.12)] * 2,
                             seed=99, max_iter=5000)
    assert r1["x"] == r2["x"]


def test_respects_bounds():
    r = simulated_annealing(lambda v: (v[0] - 100) ** 2, [0.0], bounds=[(-5, 5)],
                            max_iter=5000)
    assert abs(r["x"][0] - 5.0) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        simulated_annealing(lambda v: v[0], [])
    with pytest.raises(ValueError):
        simulated_annealing(lambda v: v[0], [1.0], cooling=1.5)
