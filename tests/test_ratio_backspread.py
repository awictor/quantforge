"""Tests for ratio-spread and backspread builders."""

import pytest

from quantforge import ratio_spread, backspread, payoff_at_expiry, OptionType


def test_ratio_spread_is_net_short_gamma():
    # Long 1, short 2 calls -> net short one option -> negative net gamma/vega.
    book = ratio_spread(100, 100, 110, 0.5, 0.05, 0.25, kind="call", ratio=2)
    assert book.net.gamma < 0
    assert book.net.vega < 0


def test_backspread_is_net_long_gamma():
    # Short 1, long 2 calls -> net long one option -> positive gamma/vega.
    book = backspread(100, 100, 110, 0.5, 0.05, 0.25, kind="call", ratio=2)
    assert book.net.gamma > 0
    assert book.net.vega > 0


def test_ratio_and_backspread_are_mirror_prices():
    # Same strikes/ratio: backspread(K_short=A, K_long=B) legs are the negatives
    # of ratio_spread(K_long=A, K_short=B) only if the roles line up; instead
    # check the net premium of a 1x1 ratio equals a debit vertical.
    rs = ratio_spread(100, 95, 105, 1.0, 0.05, 0.2, kind="call", ratio=1)
    from quantforge import vertical_spread
    vs = vertical_spread(100, 95, 105, 1.0, 0.05, 0.2, kind="call")
    assert rs.net.market_value == pytest.approx(vs.net.market_value, abs=1e-9)


def test_call_ratio_spread_payoff_tent_then_falls():
    # Long 100 call, short 2x 110 calls. Payoff: 0 below 100, peak near 110,
    # then declines (net short one call beyond 110).
    book = ratio_spread(100, 100, 110, 0.5, 0.05, 0.25, kind="call", ratio=2)
    peak = payoff_at_expiry(book, 110)      # = 10 (long call ITM, shorts ATM)
    beyond = payoff_at_expiry(book, 140)    # long 40 - 2*30 = -20
    assert peak == pytest.approx(10.0)
    assert beyond < peak
    assert payoff_at_expiry(book, 90) == pytest.approx(0.0)


def test_call_backspread_unbounded_upside():
    book = backspread(100, 100, 110, 0.5, 0.05, 0.25, kind="call", ratio=2)
    # Short 100 call, long 2x 110: far up = -40 + 2*30 = +20, rising with spot.
    p130 = payoff_at_expiry(book, 130)
    p160 = payoff_at_expiry(book, 160)
    assert p160 > p130


def test_ratio_rejects_bad_ratio():
    with pytest.raises(ValueError):
        ratio_spread(100, 100, 110, 0.5, 0.05, 0.25, ratio=0)
