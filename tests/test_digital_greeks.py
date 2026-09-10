"""Tests for cash-or-nothing digital Greeks and pin risk."""

import pytest

from quantforge import (
    digital_greeks, cash_or_nothing, spread_option, exchange_option, OptionType,
)


def test_delta_matches_re_difference():
    g = digital_greeks(100, 100, 0.5, 0.05, 0.25, OptionType.CALL)
    h = 0.05
    up = cash_or_nothing(100 + h, 100, 0.5, 0.05, 0.25, OptionType.CALL)
    dn = cash_or_nothing(100 - h, 100, 0.5, 0.05, 0.25, OptionType.CALL)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-4)


def test_call_delta_positive_put_negative():
    c = digital_greeks(100, 100, 0.5, 0.05, 0.25, OptionType.CALL)
    p = digital_greeks(100, 100, 0.5, 0.05, 0.25, OptionType.PUT)
    assert c["delta"] > 0 and p["delta"] < 0


def test_pin_risk_delta_grows_as_expiry_nears():
    # At the strike, the digital delta spikes as time to expiry shrinks.
    d_far = digital_greeks(100, 100, 0.5, 0.05, 0.25, OptionType.CALL)["delta"]
    d_near = digital_greeks(100, 100, 0.02, 0.05, 0.25, OptionType.CALL)["delta"]
    assert d_near > d_far * 2


def test_price_field_matches_direct():
    g = digital_greeks(100, 105, 0.5, 0.04, 0.3, OptionType.CALL, cash=10.0)
    assert g["price"] == pytest.approx(
        cash_or_nothing(100, 105, 0.5, 0.04, 0.3, OptionType.CALL, cash=10.0))


# --- Regression guard: Kirk spread at K=0 vs exact Margrabe ---
def test_kirk_zero_strike_matches_margrabe():
    # A spread option with strike 0 is an exchange option; Kirk should match the
    # exact Margrabe formula there (both use r-adjusted forwards internally).
    S1, S2, t, r, s1, s2, rho = 100, 95, 1.0, 0.0, 0.2, 0.25, 0.3
    kirk = spread_option(S1, S2, 0.0, t, r, s1, s2, rho, option_type=OptionType.CALL)
    margrabe = exchange_option(S1, S2, t, s1, s2, rho)
    assert kirk == pytest.approx(margrabe, abs=1e-6)
