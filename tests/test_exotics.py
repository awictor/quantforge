"""Tests for exotic option closed forms.

Strategy: exotics have few published reference tables, so we lean on exact
analytic identities that must hold regardless of parameters:
  * digital decomposition: asset_or_nothing - K*cash_or_nothing = vanilla,
  * in-out parity: knock-in + knock-out = vanilla (zero rebate),
  * knock-out below/above a far barrier -> vanilla,
  * geometric Asian < vanilla (averaging cuts volatility),
plus one Haug reference value for a barrier.
"""

import math

import pytest

from quantforge import (
    call_price, put_price, price, OptionType,
    cash_or_nothing, asset_or_nothing, barrier_option, geometric_asian, Barrier,
)


# --- Digital decomposition: C_vanilla = asset_on - K * cash_on(cash=1) ---
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
@pytest.mark.parametrize("S,K", [(100, 90), (100, 100), (100, 110)])
def test_digital_decomposition(ot, S, K):
    t, r, sigma = 1.0, 0.05, 0.25
    aon = asset_or_nothing(S, K, t, r, sigma, ot)
    con = cash_or_nothing(S, K, t, r, sigma, ot, cash=1.0)
    vanilla = price(S, K, t, r, sigma, ot)
    if ot is OptionType.CALL:
        assert aon - K * con == pytest.approx(vanilla, abs=1e-10)
    else:
        # Put: K*cash_on - asset_on = vanilla put.
        assert K * con - aon == pytest.approx(vanilla, abs=1e-10)


def test_cash_or_nothing_bounded_by_discounted_cash():
    v = cash_or_nothing(100, 100, 1.0, 0.05, 0.2, OptionType.CALL, cash=10.0)
    assert 0 < v < 10.0 * math.exp(-0.05)


# --- In-out parity: in + out = vanilla (rebate 0) ---
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
@pytest.mark.parametrize("pair", [
    (Barrier.DOWN_IN, Barrier.DOWN_OUT, 90.0),   # barrier below spot
    (Barrier.UP_IN, Barrier.UP_OUT, 120.0),      # barrier above spot
])
@pytest.mark.parametrize("K", [95, 100, 105])
def test_in_out_parity(ot, pair, K):
    b_in, b_out, H = pair
    S, t, r, sigma = 100.0, 1.0, 0.05, 0.25
    ki = barrier_option(S, K, H, t, r, sigma, ot, b_in)
    ko = barrier_option(S, K, H, t, r, sigma, ot, b_out)
    vanilla = price(S, K, t, r, sigma, ot)
    assert ki + ko == pytest.approx(vanilla, abs=1e-8)


# --- A far, un-hittable knock-out equals the vanilla option ---
def test_far_down_out_call_equals_vanilla():
    # Barrier far below spot -> practically never knocked out.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.25
    ko = barrier_option(S, K, 1.0, t, r, sigma, OptionType.CALL, Barrier.DOWN_OUT)
    assert ko == pytest.approx(call_price(S, K, t, r, sigma), abs=1e-4)


def test_far_up_out_put_equals_vanilla():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.25
    ko = barrier_option(S, K, 1e6, t, r, sigma, OptionType.PUT, Barrier.UP_OUT)
    assert ko == pytest.approx(put_price(S, K, t, r, sigma), abs=1e-4)


# --- Knock-out prices are non-negative and below vanilla ---
@pytest.mark.parametrize("bar,H", [
    (Barrier.DOWN_OUT, 90), (Barrier.UP_OUT, 120),
    (Barrier.DOWN_IN, 90), (Barrier.UP_IN, 120),
])
def test_barrier_prices_nonnegative(bar, H):
    v = barrier_option(100, 100, H, 1.0, 0.05, 0.25, OptionType.CALL, bar)
    assert v >= -1e-9


# --- Haug reference: down-and-out call with a cash rebate ---
# Haug (2007) "Complete Guide to Option Pricing Formulas": S=100, K=90, H=95,
# t=0.5, r=0.08, b=0.04, sigma=0.25, rebate=3 -> 9.0246.
def test_haug_down_out_call_reference():
    v = barrier_option(100, 90, 95, 0.5, 0.08, 0.25, OptionType.CALL,
                       Barrier.DOWN_OUT, b=0.04, rebate=3.0)
    assert v == pytest.approx(9.0246, abs=5e-3)


# --- Geometric Asian ---
def test_geometric_asian_cheaper_than_vanilla_call():
    # Averaging reduces effective vol, so the Asian call is worth less.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.3
    asian = geometric_asian(S, K, t, r, sigma, OptionType.CALL)
    vanilla = call_price(S, K, t, r, sigma)
    assert 0 < asian < vanilla


def test_geometric_asian_put_call_parity():
    # Parity holds with the Asian's adjusted carry/vol.
    S, K, t, r, sigma = 100, 95, 1.0, 0.04, 0.25
    c = geometric_asian(S, K, t, r, sigma, OptionType.CALL)
    p = geometric_asian(S, K, t, r, sigma, OptionType.PUT)
    sigma_a = sigma / math.sqrt(3.0)
    b_a = 0.5 * (r - sigma * sigma / 6.0)
    lhs = c - p
    rhs = S * math.exp((b_a - r) * t) - K * math.exp(-r * t)
    assert lhs == pytest.approx(rhs, abs=1e-10)
