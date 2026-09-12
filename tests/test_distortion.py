"""Distortion (spectral) risk pricing: Wang and proportional hazard."""

import pytest

from quantforge import (
    wang_premium, proportional_hazard_premium, expected_loss,
)
from quantforge.panjer import panjer_poisson, aggregate_mean


def _g():
    return panjer_poisson(3.0, [0.0, 0.4, 0.6])


def test_wang_zero_lambda_is_mean():
    g = _g()
    assert abs(wang_premium(g, 0.0) - aggregate_mean(g)) < 1e-4


def test_expected_loss_matches_mean():
    g = _g()
    assert abs(expected_loss(g) - aggregate_mean(g)) < 1e-9


def test_wang_monotone_in_lambda():
    g = _g()
    assert wang_premium(g, 1.0) > wang_premium(g, 0.5) > wang_premium(g, 0.0)
    assert wang_premium(g, -0.5) < aggregate_mean(g)


def test_ph_unit_rho_is_mean_and_loads_above():
    g = _g()
    assert abs(proportional_hazard_premium(g, 1.0) - aggregate_mean(g)) < 1e-4
    assert proportional_hazard_premium(g, 2.0) > aggregate_mean(g)
    assert (proportional_hazard_premium(g, 3.0)
            > proportional_hazard_premium(g, 2.0))


def test_validation():
    with pytest.raises(ValueError):
        proportional_hazard_premium(_g(), 0.5)
