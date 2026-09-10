"""Tests for implied dividend-curve bootstrapping."""

import math

import pytest

from quantforge import dividend_curve, call_price, put_price


def _chain(S, t, r, q, sigma, strikes):
    b = r - q
    calls = [call_price(S, K, t, r, sigma, b=b) for K in strikes]
    puts = [put_price(S, K, t, r, sigma, b=b) for K in strikes]
    return calls, puts


def test_recovers_known_dividend_term_structure():
    S, r, sigma = 100.0, 0.05, 0.25
    strikes = [80, 90, 100, 110, 120]
    div_by_t = {0.25: 0.01, 0.5: 0.02, 1.0: 0.03}   # rising dividend yield

    chain = []
    for t, q in div_by_t.items():
        calls, puts = _chain(S, t, r, q, sigma, strikes)
        chain.append((t, strikes, calls, puts))

    curve = dividend_curve(chain, spot=S)
    # Sorted by expiry.
    assert [t for t, _ in curve] == [0.25, 0.5, 1.0]
    for t, res in curve:
        assert res.implied_div_yield == pytest.approx(div_by_t[t], abs=1e-6)
        assert res.implied_rate == pytest.approx(r, abs=1e-6)


def test_forward_below_spot_for_high_dividend():
    S, r, sigma = 100.0, 0.02, 0.25
    strikes = [80, 90, 100, 110, 120]
    calls, puts = _chain(S, 1.0, r, 0.08, sigma, strikes)  # q > r
    curve = dividend_curve([(1.0, strikes, calls, puts)], spot=S)
    _, res = curve[0]
    assert res.forward < S
    assert res.implied_div_yield == pytest.approx(0.08, abs=1e-6)


def test_curve_is_sorted_by_expiry():
    S, r, sigma = 100.0, 0.05, 0.25
    strikes = [90, 100, 110]
    chain = []
    for t in (1.0, 0.25, 0.5):   # deliberately unsorted
        calls, puts = _chain(S, t, r, 0.02, sigma, strikes)
        chain.append((t, strikes, calls, puts))
    curve = dividend_curve(chain, spot=S)
    ts = [t for t, _ in curve]
    assert ts == sorted(ts)
