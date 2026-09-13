"""Differential evolution global optimization."""

import math

import pytest

from quantforge import differential_evolution, nelder_mead


def test_sphere():
    r = differential_evolution(lambda v: sum(x * x for x in v), [(-5, 5)] * 3)
    assert r["fun"] < 1e-6
    assert all(abs(x) < 1e-3 for x in r["x"])


def test_rosenbrock():
    rosen = lambda v: (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2
    r = differential_evolution(rosen, [(-5, 5), (-5, 5)], max_iter=2000)
    assert abs(r["x"][0] - 1) < 1e-2 and abs(r["x"][1] - 1) < 1e-2


def test_beats_local_method_on_multimodal():
    rastrigin = lambda v: 10 * len(v) + sum(x * x - 10 * math.cos(2 * math.pi * x) for x in v)
    de = differential_evolution(rastrigin, [(-5.12, 5.12)] * 2, max_iter=3000)
    _, nm_f = nelder_mead(rastrigin, [3.0, 3.0], step=0.5)
    assert de["fun"] < 1e-3           # DE finds the global min
    assert de["fun"] < nm_f           # local method gets stuck higher


def test_deterministic_for_seed():
    rosen = lambda v: (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2
    r1 = differential_evolution(rosen, [(-5, 5), (-5, 5)], seed=42, max_iter=500)
    r2 = differential_evolution(rosen, [(-5, 5), (-5, 5)], seed=42, max_iter=500)
    assert r1["x"] == r2["x"]


def test_respects_bounds():
    r = differential_evolution(lambda v: (v[0] - 100) ** 2, [(-5, 5)])
    assert abs(r["x"][0] - 5.0) < 1e-6      # true min at 100, clamped to upper bound


def test_ackley():
    def ackley(v):
        a, b, c, n = 20, 0.2, 2 * math.pi, len(v)
        return (-a * math.exp(-b * math.sqrt(sum(x * x for x in v) / n))
                - math.exp(sum(math.cos(c * x) for x in v) / n) + a + math.e)
    r = differential_evolution(ackley, [(-5, 5)] * 3, max_iter=3000)
    assert r["fun"] < 1e-3


def test_validation():
    with pytest.raises(ValueError):
        differential_evolution(lambda v: v[0], [])
    with pytest.raises(ValueError):
        differential_evolution(lambda v: v[0], [(5, -5)])       # lo > hi
    with pytest.raises(ValueError):
        differential_evolution(lambda v: v[0], [(-5, 5)], pop_size=2)
