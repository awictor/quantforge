"""Tests for covered-call and cash-secured-put income analytics."""

import pytest

from quantforge import covered_call, cash_secured_put, call_price, put_price


def test_covered_call_premium_is_call_value():
    m = covered_call(100, 105, 0.25, 0.04, 0.25)
    assert m.premium == pytest.approx(call_price(100, 105, 0.25, 0.04, 0.25))


def test_covered_call_breakeven_and_yield():
    m = covered_call(100, 105, 0.25, 0.04, 0.25)
    assert m.breakeven == pytest.approx(100 - m.premium)
    assert m.static_yield == pytest.approx(m.premium / 100)
    assert m.annualized_yield == pytest.approx(m.static_yield / 0.25)


def test_covered_call_if_assigned_return():
    S, K = 100, 105
    m = covered_call(S, K, 0.25, 0.04, 0.25)
    assert m.if_assigned_return == pytest.approx(((K - S) + m.premium) / S)


def test_cash_secured_put_premium_is_put_value():
    m = cash_secured_put(100, 95, 0.25, 0.04, 0.3)
    assert m.premium == pytest.approx(put_price(100, 95, 0.25, 0.04, 0.3))


def test_cash_secured_put_breakeven_is_strike_minus_premium():
    m = cash_secured_put(100, 95, 0.25, 0.04, 0.3)
    assert m.breakeven == pytest.approx(95 - m.premium)
    assert m.static_yield == pytest.approx(m.premium / 95)


def test_custom_premium_overrides_model():
    m = covered_call(100, 105, 0.25, 0.04, 0.25, premium=4.0)
    assert m.premium == 4.0
    assert m.breakeven == pytest.approx(96.0)


def test_annualized_yield_scales_with_time():
    short = covered_call(100, 105, 0.25, 0.04, 0.25)
    long = covered_call(100, 105, 1.0, 0.04, 0.25)
    # A shorter tenor annualizes a similar-ish premium to a higher yield.
    assert short.annualized_yield > long.annualized_yield
