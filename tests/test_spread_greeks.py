"""Tests for Kirk spread-option Greeks (finite differences)."""

import pytest

from quantforge import spread_greeks, spread_option, OptionType


BASE = dict(S1=100, S2=95, K=5, t=1.0, r=0.05, sigma1=0.2, sigma2=0.25, rho=0.5)


def test_delta_signs():
    g = spread_greeks(**BASE)
    assert g["delta1"] > 0
    assert g["delta2"] < 0


def test_delta1_matches_re_difference():
    g = spread_greeks(**BASE)
    h = 0.05
    up = spread_option(100 + h, 95, 5, 1.0, 0.05, 0.2, 0.25, 0.5)
    dn = spread_option(100 - h, 95, 5, 1.0, 0.05, 0.2, 0.25, 0.5)
    assert g["delta1"] == pytest.approx((up - dn) / (2 * h), abs=1e-3)


def test_cross_gamma_near_negative_own_gamma():
    g = spread_greeks(**BASE)
    assert g["cross"] == pytest.approx(-g["gamma1"], abs=2e-3)


def test_correlation_vega_negative():
    g = spread_greeks(**BASE)
    assert g["corr_vega"] < 0


def test_own_gammas_positive():
    g = spread_greeks(**BASE)
    assert g["gamma1"] > 0 and g["gamma2"] > 0


def test_price_field_matches_direct():
    g = spread_greeks(**BASE)
    assert g["price"] == pytest.approx(
        spread_option(100, 95, 5, 1.0, 0.05, 0.2, 0.25, 0.5, option_type=OptionType.CALL),
        abs=1e-9)
