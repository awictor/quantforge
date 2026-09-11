"""Implied discount factor from two-strike call-put pairs (bsm module)."""

import math

import pytest

from quantforge import (
    implied_discount_factor, implied_forward_from_parity,
    call_price, put_price, forward_price,
)


S, T, R, SIG = 100.0, 1.0, 0.05, 0.2
K1, K2 = 90.0, 110.0


def _pairs(k, b=None):
    return call_price(S, k, T, R, SIG, b=b), put_price(S, k, T, R, SIG, b=b)


def test_recovers_discount_factor():
    c1, p1 = _pairs(K1)
    c2, p2 = _pairs(K2)
    df = implied_discount_factor(c1, p1, K1, c2, p2, K2)
    assert df == pytest.approx(math.exp(-R * T), abs=1e-9)


def test_recovers_with_dividend():
    b = R - 0.03
    c1, p1 = _pairs(K1, b=b)
    c2, p2 = _pairs(K2, b=b)
    # Discounting still uses r regardless of carry.
    df = implied_discount_factor(c1, p1, K1, c2, p2, K2)
    assert df == pytest.approx(math.exp(-R * T), abs=1e-9)


def test_strike_order_independent():
    c1, p1 = _pairs(K1)
    c2, p2 = _pairs(K2)
    a = implied_discount_factor(c1, p1, K1, c2, p2, K2)
    b = implied_discount_factor(c2, p2, K2, c1, p1, K1)
    assert a == pytest.approx(b, abs=1e-12)


def test_combined_with_implied_forward():
    # DF and forward together reconstruct each pair's C - P.
    c1, p1 = _pairs(K1)
    c2, p2 = _pairs(K2)
    df = implied_discount_factor(c1, p1, K1, c2, p2, K2)
    F = implied_forward_from_parity(c1, p1, K1, T, R)
    assert F == pytest.approx(forward_price(S, T, R), abs=1e-9)
    assert df * (F - K2) == pytest.approx(c2 - p2, abs=1e-9)


def test_same_strike_raises():
    c1, p1 = _pairs(K1)
    with pytest.raises(ValueError):
        implied_discount_factor(c1, p1, K1, c1, p1, K1)
