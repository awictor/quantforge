"""Tests for implied-forward / dividend extraction from an option chain."""

import math

import pytest

from quantforge import implied_forward, call_price, put_price


def _synthetic_chain(S, t, r, q, sigma, strikes):
    """Build a European call/put chain from BSM with carry b = r - q."""
    b = r - q
    calls = [call_price(S, K, t, r, sigma, b=b) for K in strikes]
    puts = [put_price(S, K, t, r, sigma, b=b) for K in strikes]
    return calls, puts


def test_recovers_forward_and_discount():
    S, t, r, q, sigma = 100.0, 1.0, 0.05, 0.02, 0.25
    strikes = [80, 90, 100, 110, 120]
    calls, puts = _synthetic_chain(S, t, r, q, sigma, strikes)

    res = implied_forward(strikes, calls, puts, t, spot=S)

    F_true = S * math.exp((r - q) * t)
    D_true = math.exp(-r * t)
    assert res.forward == pytest.approx(F_true, abs=1e-6)
    assert res.discount_factor == pytest.approx(D_true, abs=1e-8)
    assert res.implied_rate == pytest.approx(r, abs=1e-6)
    assert res.implied_div_yield == pytest.approx(q, abs=1e-6)
    assert res.rmse < 1e-6


def test_no_dividend_forward():
    S, t, r, q, sigma = 50.0, 0.5, 0.03, 0.0, 0.4
    strikes = [40, 45, 50, 55, 60]
    calls, puts = _synthetic_chain(S, t, r, q, sigma, strikes)
    res = implied_forward(strikes, calls, puts, t, spot=S)
    assert res.forward == pytest.approx(S * math.exp(r * t), abs=1e-6)
    assert res.implied_div_yield == pytest.approx(0.0, abs=1e-6)


def test_high_dividend_forward_below_spot():
    # Dividend yield above the rate -> forward below spot.
    S, t, r, q, sigma = 100.0, 1.0, 0.03, 0.08, 0.2
    strikes = [80, 90, 100, 110, 120]
    calls, puts = _synthetic_chain(S, t, r, q, sigma, strikes)
    res = implied_forward(strikes, calls, puts, t, spot=S)
    assert res.forward < S
    assert res.implied_div_yield == pytest.approx(q, abs=1e-6)


def test_requires_two_strikes():
    with pytest.raises(ValueError):
        implied_forward([100], [5.0], [4.0], t=1.0)


def test_mismatched_lengths_rejected():
    with pytest.raises(ValueError):
        implied_forward([90, 100], [5.0], [4.0, 3.0], t=1.0)


def test_robust_to_small_noise():
    S, t, r, q, sigma = 100.0, 1.0, 0.05, 0.02, 0.25
    strikes = [70, 85, 100, 115, 130]
    calls, puts = _synthetic_chain(S, t, r, q, sigma, strikes)
    # Add tiny symmetric bid/ask noise; least-squares should stay close.
    noise = [0.01, -0.01, 0.005, -0.005, 0.0]
    calls = [c + n for c, n in zip(calls, noise)]
    puts = [p - n for p, n in zip(puts, noise)]
    res = implied_forward(strikes, calls, puts, t, spot=S)
    F_true = S * math.exp((r - q) * t)
    assert res.forward == pytest.approx(F_true, rel=0.02)
