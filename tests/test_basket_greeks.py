"""Tests for two-asset basket-option Greeks (finite differences)."""

import pytest

from quantforge import basket_greeks, basket_option, OptionType


BASE = dict(spots=(100, 100), weights=(0.5, 0.5), K=100, t=1.0, r=0.05,
            sigmas=(0.2, 0.3), corr=0.4)


def test_both_deltas_positive_for_call():
    g = basket_greeks(**BASE)
    assert g["delta1"] > 0 and g["delta2"] > 0


def test_delta1_matches_re_difference():
    g = basket_greeks(**BASE)
    h = 0.05
    up = basket_option((100 + h, 100), (0.5, 0.5), 100, 1.0, 0.05, (0.2, 0.3), 0.4)
    dn = basket_option((100 - h, 100), (0.5, 0.5), 100, 1.0, 0.05, (0.2, 0.3), 0.4)
    assert g["delta1"] == pytest.approx((up - dn) / (2 * h), abs=1e-3)


def test_correlation_vega_positive():
    # Higher correlation raises the basket vol, richening the call.
    g = basket_greeks(**BASE)
    assert g["corr_vega"] > 0


def test_own_gammas_positive():
    g = basket_greeks(**BASE)
    assert g["gamma1"] > 0 and g["gamma2"] > 0


def test_price_field_matches_direct():
    g = basket_greeks(**BASE)
    assert g["price"] == pytest.approx(
        basket_option((100, 100), (0.5, 0.5), 100, 1.0, 0.05, (0.2, 0.3), 0.4,
                      option_type=OptionType.CALL), abs=1e-9)


def test_weight_shifts_delta_balance():
    # Overweighting asset 1 raises its delta relative to asset 2.
    g = basket_greeks(spots=(100, 100), weights=(0.8, 0.2), K=100, t=1.0, r=0.05,
                      sigmas=(0.2, 0.3), corr=0.4)
    assert g["delta1"] > g["delta2"]
