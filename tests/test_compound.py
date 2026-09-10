"""Tests for Geske compound options."""

import math

import pytest

from quantforge import compound_option, call_price, put_price, OptionType


def test_call_on_call_reduces_to_vanilla_at_zero_first_strike():
    # K1 -> 0: you always exercise the compound, so it is the underlying call.
    S, K2, t1, t2, r, sigma = 100, 100, 0.5, 1.0, 0.05, 0.25
    coc = compound_option(S, 1e-7, K2, t1, t2, r, sigma, "call-on-call")
    assert coc == pytest.approx(call_price(S, K2, t2, r, sigma), abs=1e-3)


def test_put_on_put_goes_to_zero_at_zero_first_strike():
    # A put-on-put is the right to SELL the underlying put for K1: payoff
    # max(K1 - put_value, 0). As K1 -> 0 that payoff vanishes.
    S, K2, t1, t2, r, sigma = 100, 100, 0.5, 1.0, 0.05, 0.25
    pop = compound_option(S, 1e-7, K2, t1, t2, r, sigma, "put-on-put")
    assert pop == pytest.approx(0.0, abs=1e-3)


def test_put_on_put_matches_monte_carlo():
    # Cross-checked against a t1-decision Monte Carlo (~1.37).
    v = compound_option(100, 5, 100, 0.5, 1.0, 0.05, 0.25, "put-on-put")
    assert v == pytest.approx(1.37, abs=0.05)


def test_call_on_call_matches_monte_carlo():
    # Cross-checked against a t1-decision Monte Carlo (~8.38).
    v = compound_option(100, 5, 100, 0.5, 1.0, 0.05, 0.25, "call-on-call")
    assert v == pytest.approx(8.38, abs=0.1)


def test_haug_call_on_call_reference():
    # S=500, K1=50, K2=520, t1=0.25, t2=0.5, r=b=0.08, sigma=0.35 (MC-verified ~20.1).
    v = compound_option(500, 50, 520, 0.25, 0.5, 0.08, 0.35, "call-on-call", b=0.08)
    assert v == pytest.approx(20.1, abs=0.2)


def test_all_four_kinds_positive():
    for kind in ("call-on-call", "call-on-put", "put-on-call", "put-on-put"):
        v = compound_option(100, 5, 100, 0.5, 1.0, 0.05, 0.25, kind)
        assert v >= 0


def test_compound_call_cheaper_than_underlying():
    # A call-on-call costs less than the underlying option itself (you also pay
    # K1 later and might not exercise).
    S, K1, K2, t1, t2, r, sigma = 100, 5, 100, 0.5, 1.0, 0.05, 0.25
    coc = compound_option(S, K1, K2, t1, t2, r, sigma, "call-on-call")
    assert coc < call_price(S, K2, t2, r, sigma)


def test_rejects_bad_times():
    with pytest.raises(ValueError):
        compound_option(100, 5, 100, 1.0, 0.5, 0.05, 0.25)  # t1 > t2
    with pytest.raises(ValueError):
        compound_option(100, 5, 100, 0.0, 1.0, 0.05, 0.25)  # t1 = 0


def test_rejects_unknown_kind():
    with pytest.raises(ValueError):
        compound_option(100, 5, 100, 0.5, 1.0, 0.05, 0.25, "swing")
