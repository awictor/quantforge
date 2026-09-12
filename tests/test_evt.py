"""Extreme value theory: Hill index and peaks-over-threshold VaR / ES."""

import math
import random

import pytest

from quantforge import (
    hill_estimator, gpd_fit_pot, gpd_var, gpd_expected_shortfall,
    gev_cdf, gev_return_level, gev_fit_block_maxima,
)


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


def test_gev_cdf_monotone_and_bounded():
    xs = [-2, -1, 0, 1, 2, 3, 5]
    cdfs = [gev_cdf(x, 0, 1, 0.2) for x in xs]
    assert all(cdfs[i] <= cdfs[i + 1] for i in range(len(cdfs) - 1))
    assert all(0 <= c <= 1 for c in cdfs)


def test_gev_gumbel_limit():
    assert gev_cdf(1.0, 0, 1, 1e-10) == pytest.approx(math.exp(-math.exp(-1.0)), abs=1e-6)


def test_return_level_increases_with_period():
    assert gev_return_level(100, 0, 1, 0.2) > gev_return_level(10, 0, 1, 0.2) \
        > gev_return_level(2, 0, 1, 0.2)


def test_return_level_is_quantile():
    T = 50
    lev = gev_return_level(T, 0, 1, 0.2)
    assert gev_cdf(lev, 0, 1, 0.2) == pytest.approx(1 - 1 / T, abs=1e-9)


def test_gumbel_fit_recovers_parameters():
    random.seed(3)
    bm = [-math.log(-math.log(random.random())) for _ in range(5000)]
    loc, scale, shape = gev_fit_block_maxima(bm)
    assert abs(loc) < 0.1
    assert abs(scale - 1) < 0.1
    assert shape == 0.0


def test_gev_validation():
    with pytest.raises(ValueError):
        gev_return_level(1, 0, 1, 0.2)
    with pytest.raises(ValueError):
        gev_cdf(1, 0, 0, 0.2)


def test_validation():
    losses = _pareto_losses(3.0, 5000)
    with pytest.raises(ValueError):
        hill_estimator(losses, 0)
    with pytest.raises(ValueError):
        gpd_fit_pot(losses, 1e9)   # no exceedances
