"""Tests for option-strategy builders."""

import math

import pytest

from quantforge import (
    vertical_spread, straddle, strangle, risk_reversal, butterfly, iron_condor,
    payoff_at_expiry, payoff_profile, break_evens, price_book, Contract,
    call_price, put_price, OptionType,
)


BASE = dict(t=0.5, r=0.04, sigma=0.25)


def test_bull_call_spread_payoff_caps():
    # Long 95 call, short 105 call: payoff 0 below 95, ramps to 10 above 105.
    book = vertical_spread(100, 95, 105, kind="call", **BASE)
    assert payoff_at_expiry(book, 90) == pytest.approx(0.0)
    assert payoff_at_expiry(book, 100) == pytest.approx(5.0)
    assert payoff_at_expiry(book, 120) == pytest.approx(10.0)  # capped


def test_bull_call_spread_net_premium_is_debit():
    # Long the lower strike (more expensive) => net premium paid > 0.
    book = vertical_spread(100, 95, 105, kind="call", **BASE)
    assert book.net.market_value > 0


def test_straddle_payoff_is_v_shaped():
    book = straddle(100, 100, **BASE)
    assert payoff_at_expiry(book, 100) == pytest.approx(0.0)
    assert payoff_at_expiry(book, 120) == pytest.approx(20.0)
    assert payoff_at_expiry(book, 80) == pytest.approx(20.0)


def test_straddle_break_evens_straddle_the_strike():
    book = straddle(100, 100, **BASE)
    prem = book.net.market_value
    bes = break_evens(book, 50, 150)
    assert len(bes) == 2
    lo, hi = sorted(bes)
    # Break-evens are strike +/- premium.
    assert lo == pytest.approx(100 - prem, abs=0.1)
    assert hi == pytest.approx(100 + prem, abs=0.1)


def test_strangle_flat_between_strikes():
    book = strangle(100, 90, 110, **BASE)
    assert payoff_at_expiry(book, 100) == pytest.approx(0.0)  # between strikes
    assert payoff_at_expiry(book, 85) == pytest.approx(5.0)   # put ITM
    assert payoff_at_expiry(book, 115) == pytest.approx(5.0)  # call ITM


def test_risk_reversal_is_long_delta():
    # Short put + long call => positive net delta (synthetic long-ish).
    book = risk_reversal(100, 90, 110, **BASE)
    assert book.net.delta > 0


def test_butterfly_payoff_peaks_at_middle_strike():
    book = butterfly(100, 90, 100, 110, kind="call", **BASE)
    # Peak payoff at the middle strike = wing width.
    assert payoff_at_expiry(book, 100) == pytest.approx(10.0)
    assert payoff_at_expiry(book, 90) == pytest.approx(0.0)
    assert payoff_at_expiry(book, 110) == pytest.approx(0.0)
    assert payoff_at_expiry(book, 130) == pytest.approx(0.0)  # wings cancel


def test_butterfly_is_cheap_debit():
    book = butterfly(100, 90, 100, 110, kind="call", **BASE)
    # Net premium is a small positive debit, below the max payoff.
    assert 0 < book.net.market_value < 10


def test_iron_condor_collects_premium():
    # Selling both spreads => net credit (negative market value to us).
    book = iron_condor(100, 80, 90, 110, 120, **BASE)
    assert book.net.market_value < 0
    # Max profit region: flat between the short strikes.
    assert payoff_at_expiry(book, 100) == pytest.approx(0.0)
    # Wings cap the loss.
    assert payoff_at_expiry(book, 70) == pytest.approx(-10.0)
    assert payoff_at_expiry(book, 130) == pytest.approx(-10.0)


def test_payoff_profile_length():
    book = straddle(100, 100, **BASE)
    spots = [80, 90, 100, 110, 120]
    prof = payoff_profile(book, spots)
    assert len(prof) == 5


def test_net_price_matches_leg_sum():
    book = vertical_spread(100, 95, 105, kind="call", **BASE)
    long_c = call_price(100, 95, 0.5, 0.04, 0.25)
    short_c = call_price(100, 105, 0.5, 0.04, 0.25)
    assert book.net.market_value == pytest.approx(long_c - short_c, abs=1e-9)
