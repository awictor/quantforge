"""Tests for the strategy P&L summary report."""

import pytest

from quantforge import (
    strategy_report, iron_condor, straddle, vertical_spread, backspread,
)


def test_iron_condor_bounded_both_ways():
    r = strategy_report(iron_condor(100, 80, 90, 110, 120, 0.5, 0.04, 0.25))
    assert not r["profit_unbounded"]
    assert not r["loss_unbounded"]
    assert len(r["break_evens"]) == 2
    assert r["max_profit"] > 0 > r["max_loss"]


def test_straddle_has_unbounded_profit_capped_loss():
    r = strategy_report(straddle(100, 100, 0.5, 0.04, 0.25))
    assert r["profit_unbounded"]
    assert not r["loss_unbounded"]
    assert len(r["break_evens"]) == 2


def test_vertical_spread_caps_both():
    r = strategy_report(vertical_spread(100, 95, 105, 0.5, 0.04, 0.25, kind="call"))
    assert not r["profit_unbounded"]
    assert not r["loss_unbounded"]
    # Max profit + max loss ~ the strike width (10) minus/plus premium symmetry.
    assert r["max_profit"] + abs(r["max_loss"]) == pytest.approx(10.0, abs=0.05)


def test_call_backspread_unbounded_profit():
    r = strategy_report(backspread(100, 100, 110, 0.5, 0.05, 0.25, kind="call", ratio=2))
    assert r["profit_unbounded"]


def test_net_premium_matches_book():
    from quantforge import price_book, Contract, OptionType
    book = straddle(100, 100, 0.5, 0.04, 0.25)
    r = strategy_report(book)
    assert r["net_premium"] == pytest.approx(book.net.market_value)


def test_max_loss_at_least_negative_premium_for_long_debit():
    # A long straddle's worst case is losing the full premium (spot pinned at K).
    book = straddle(100, 100, 0.5, 0.04, 0.25)
    r = strategy_report(book)
    assert r["max_loss"] == pytest.approx(-book.net.market_value, abs=0.05)
