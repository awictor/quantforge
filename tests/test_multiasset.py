"""Tests for two-asset options: exchange, spread, and basket."""

import math

import pytest

from quantforge import (
    exchange_option, spread_option, basket_option, call_price, OptionType,
)


# --- Margrabe exchange ---
def test_exchange_matches_monte_carlo():
    # Cross-checked against a correlated-GBM Monte Carlo (CF 13.13, MC ~13.13).
    v = exchange_option(100, 95, 1.0, 0.2, 0.25, rho=0.3)
    assert v == pytest.approx(13.13, abs=0.05)


def test_exchange_zero_expiry_is_intrinsic():
    assert exchange_option(110, 100, 0.0, 0.2, 0.25, 0.3) == pytest.approx(10.0)
    assert exchange_option(90, 100, 0.0, 0.2, 0.25, 0.3) == pytest.approx(0.0)


def test_exchange_higher_with_lower_correlation():
    # Lower correlation -> higher spread volatility -> more valuable exchange.
    hi_corr = exchange_option(100, 100, 1.0, 0.2, 0.2, rho=0.9)
    lo_corr = exchange_option(100, 100, 1.0, 0.2, 0.2, rho=-0.5)
    assert lo_corr > hi_corr


def test_exchange_reduces_to_zero_strike_call_when_s2_certain():
    # If asset 2 has zero vol and zero yield it is a fixed strike S2: the
    # exchange becomes a call on S1 struck at S2 (with sigma1 only).
    S1, S2, t, sigma1 = 100, 90, 1.0, 0.25
    v = exchange_option(S1, S2, t, sigma1, 0.0, rho=0.0)
    # Margrabe is rate-independent; compare to a BSM call with r=0, K=S2.
    assert v == pytest.approx(call_price(S1, S2, t, 0.0, sigma1), abs=1e-9)


# --- Kirk spread ---
def test_spread_matches_monte_carlo():
    v = spread_option(100, 95, 5, 1.0, 0.05, 0.2, 0.25, rho=0.5)
    assert v == pytest.approx(8.92, abs=0.1)


def test_spread_put_call_parity():
    kw = dict(S1=100, S2=95, K=5, t=1.0, r=0.05, sigma1=0.2, sigma2=0.25, rho=0.5)
    c = spread_option(**kw, option_type=OptionType.CALL)
    p = spread_option(**kw, option_type=OptionType.PUT)
    disc = math.exp(-0.05 * 1.0)
    F1 = 100 * math.exp(0.05); F2 = 95 * math.exp(0.05)
    assert c - p == pytest.approx(disc * (F1 - F2 - 5), abs=1e-9)


def test_spread_zero_expiry_intrinsic():
    assert spread_option(100, 90, 5, 0.0, 0.05, 0.2, 0.25, 0.3) == pytest.approx(5.0)


# --- Basket ---
def test_basket_single_asset_equals_bsm():
    # Zero weight on the second asset -> plain BSM call on the first.
    v = basket_option((100, 1), (1.0, 0.0), 100, 1.0, 0.05, (0.2, 0.2), 0.0)
    assert v == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2), abs=1e-6)


def test_basket_matches_monte_carlo():
    v = basket_option((100, 100), (0.5, 0.5), 100, 1.0, 0.05, (0.2, 0.3), 0.4)
    assert v == pytest.approx(10.85, abs=0.1)


def test_basket_put_call_parity():
    kw = dict(spots=(100, 100), weights=(0.5, 0.5), K=100, t=1.0, r=0.05,
              sigmas=(0.2, 0.3), corr=0.4)
    c = basket_option(**kw, option_type=OptionType.CALL)
    p = basket_option(**kw, option_type=OptionType.PUT)
    # Forward of the basket = 0.5*100*e^{rt} + 0.5*100*e^{rt} = 100 e^{rt}.
    disc = math.exp(-0.05)
    M1 = 100 * math.exp(0.05)
    assert c - p == pytest.approx(disc * (M1 - 100), abs=1e-9)


def test_basket_rejects_wrong_asset_count():
    with pytest.raises(ValueError):
        basket_option((100, 100, 100), (0.3, 0.3, 0.4), 100, 1.0, 0.05,
                      (0.2, 0.2, 0.2), 0.0)
