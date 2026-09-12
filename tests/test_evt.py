"""Extreme value theory: Hill index and peaks-over-threshold VaR / ES."""

import random

import pytest

from quantforge import hill_estimator, gpd_fit_pot, gpd_var, gpd_expected_shortfall


def _pareto_losses(alpha, n, seed=11):
    random.seed(seed)
    return [(1.0 / (1.0 - random.random())) ** (1.0 / alpha) for _ in range(n)]


def test_hill_recovers_tail_index():
    losses = _pareto_losses(3.0, 20000)   # xi = 1/alpha = 1/3
    assert hill_estimator(losses, 1000) == pytest.approx(1 / 3, abs=0.05)


def test_gpd_positive_shape_for_heavy_tail():
    losses = _pareto_losses(3.0, 20000)
    thr = sorted(losses)[int(0.95 * len(losses))]
    xi, beta, nu, n = gpd_fit_pot(losses, thr)
    assert xi > 0


def test_es_at_least_var():
    losses = _pareto_losses(3.0, 20000)
    thr = sorted(losses)[int(0.95 * len(losses))]
    assert gpd_expected_shortfall(losses, thr, 0.99) > gpd_var(losses, thr, 0.99)


def test_var_increases_with_confidence():
    losses = _pareto_losses(3.0, 20000)
    thr = sorted(losses)[int(0.95 * len(losses))]
    assert gpd_var(losses, thr, 0.999) > gpd_var(losses, thr, 0.99)


def test_validation():
    losses = _pareto_losses(3.0, 5000)
    with pytest.raises(ValueError):
        hill_estimator(losses, 0)
    with pytest.raises(ValueError):
        gpd_fit_pot(losses, 1e9)   # no exceedances
