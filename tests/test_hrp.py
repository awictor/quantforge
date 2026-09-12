"""Hierarchical Risk Parity and inverse-volatility allocation."""

import pytest

from quantforge import inverse_volatility_weights, hierarchical_risk_parity


def test_inverse_vol_sums_and_ranking():
    cov = [[0.04, 0.01, 0.005], [0.01, 0.09, 0.01], [0.005, 0.01, 0.0025]]
    w = inverse_volatility_weights(cov)
    assert sum(w) == pytest.approx(1.0)
    assert w[2] > w[0] > w[1]   # lowest vol highest weight


def test_hrp_sums_and_positive():
    cov = [[0.04, 0.01, 0.005], [0.01, 0.09, 0.01], [0.005, 0.01, 0.0025]]
    h = hierarchical_risk_parity(cov)
    assert sum(h) == pytest.approx(1.0)
    assert all(x > 0 for x in h)


def test_hrp_block_correlated():
    # Assets 0,1 correlated (rho 0.9), asset 2 independent.
    cov = [[0.04, 0.036, 0.0], [0.036, 0.04, 0.0], [0.0, 0.0, 0.04]]
    h = hierarchical_risk_parity(cov)
    assert sum(h) == pytest.approx(1.0)
    assert h[2] > h[0] and h[2] > h[1]   # independent asset gets more


def test_hrp_equal_uncorrelated():
    cov = [[0.04, 0, 0, 0], [0, 0.04, 0, 0], [0, 0, 0.04, 0], [0, 0, 0, 0.04]]
    h = hierarchical_risk_parity(cov)
    assert all(x == pytest.approx(0.25, abs=0.05) for x in h)


def test_validation():
    with pytest.raises(ValueError):
        inverse_volatility_weights([[0, 0], [0, 0.04]])
    with pytest.raises(ValueError):
        hierarchical_risk_parity([[0, 0], [0, 0.04]])
