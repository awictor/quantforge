"""Parametric portfolio VaR / CVaR / component VaR (portopt module)."""

import math

import pytest

from quantforge import (
    portfolio_var, portfolio_cvar, component_var, portfolio_variance,
)


COV = [[0.04, 0.01, 0.0], [0.01, 0.09, 0.02], [0.0, 0.02, 0.16]]
W = [0.4, 0.4, 0.2]
Z95 = 1.6448536269514722


def test_var_matches_z_times_sigma():
    sd = math.sqrt(portfolio_variance(W, COV))
    assert portfolio_var(W, COV) == pytest.approx(Z95 * sd, abs=1e-9)


def test_cvar_exceeds_var():
    assert portfolio_cvar(W, COV) > portfolio_var(W, COV)


def test_positive_mean_lowers_var():
    with_mu = portfolio_var(W, COV, [0.05, 0.05, 0.05])
    assert with_mu < portfolio_var(W, COV)


def test_component_var_sums_to_sigma():
    sd = math.sqrt(portfolio_variance(W, COV))
    comps = component_var(W, COV)
    assert sum(comps) == pytest.approx(sd, abs=1e-9)


def test_horizon_sqrt_scaling():
    # No-mean VaR scales with sqrt(horizon): 4x horizon -> 2x VaR.
    base = portfolio_var(W, COV)
    quad = portfolio_var(W, COV, None, 0.95, 4.0)
    assert quad == pytest.approx(2.0 * base, abs=1e-9)


def test_higher_confidence_larger_var():
    assert portfolio_var(W, COV, None, 0.99) > portfolio_var(W, COV, None, 0.95)


def test_validation():
    with pytest.raises(ValueError):
        component_var([0.0, 0.0, 0.0], COV)   # zero-variance portfolio
