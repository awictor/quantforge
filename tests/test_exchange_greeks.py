"""Tests for Margrabe exchange-option Greeks (finite differences)."""

import pytest

from quantforge import exchange_greeks, exchange_option


def test_delta_signs():
    g = exchange_greeks(100, 100, 1.0, 0.2, 0.25, 0.3)
    assert g["delta1"] > 0     # long asset 1
    assert g["delta2"] < 0     # short asset 2


def test_delta1_matches_re_difference():
    g = exchange_greeks(100, 100, 1.0, 0.2, 0.25, 0.3)
    h = 0.05
    up = exchange_option(100 + h, 100, 1.0, 0.2, 0.25, 0.3)
    dn = exchange_option(100 - h, 100, 1.0, 0.2, 0.25, 0.3)
    assert g["delta1"] == pytest.approx((up - dn) / (2 * h), abs=1e-3)


def test_cross_gamma_is_negative_own_gamma():
    # Margrabe is homogeneous of degree 1, so cross-gamma = -gamma1 = -gamma2.
    g = exchange_greeks(100, 100, 1.0, 0.2, 0.25, 0.3)
    assert g["cross"] == pytest.approx(-g["gamma1"], abs=1e-4)
    assert g["gamma1"] == pytest.approx(g["gamma2"], abs=1e-4)


def test_correlation_vega_negative():
    # Higher correlation lowers the spread vol -> cheaper exchange option.
    g = exchange_greeks(100, 100, 1.0, 0.2, 0.25, 0.3)
    assert g["corr_vega"] < 0


def test_price_field_matches_direct():
    g = exchange_greeks(120, 100, 0.5, 0.3, 0.2, -0.2, q1=0.01, q2=0.02)
    assert g["price"] == pytest.approx(
        exchange_option(120, 100, 0.5, 0.3, 0.2, -0.2, 0.01, 0.02), abs=1e-9)


def test_own_gammas_positive():
    g = exchange_greeks(100, 100, 1.0, 0.2, 0.25, 0.3)
    assert g["gamma1"] > 0
    assert g["gamma2"] > 0
