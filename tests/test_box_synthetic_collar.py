"""Box spread, synthetic forward, and collar builders (strategy module)."""

import math

import pytest

from quantforge import box_spread, synthetic_forward, collar
from quantforge.strategy import payoff_at_expiry
from quantforge.bsm import price


S, R, SIG, T = 100.0, 0.05, 0.2, 1.0


def test_box_is_discounted_strike_width():
    b = box_spread(S, 90.0, 110.0, T, R, SIG)
    assert b.net.market_value == pytest.approx(math.exp(-R * T) * 20.0, abs=1e-9)


def test_box_payoff_constant():
    b = box_spread(S, 90.0, 110.0, T, R, SIG)
    assert payoff_at_expiry(b, 80.0) == pytest.approx(20.0, abs=1e-9)
    assert payoff_at_expiry(b, 100.0) == pytest.approx(20.0, abs=1e-9)
    assert payoff_at_expiry(b, 130.0) == pytest.approx(20.0, abs=1e-9)


def test_box_greeks_are_flat():
    b = box_spread(S, 90.0, 110.0, T, R, SIG)
    assert b.net.delta == pytest.approx(0.0, abs=1e-6)
    assert b.net.gamma == pytest.approx(0.0, abs=1e-8)
    assert b.net.vega == pytest.approx(0.0, abs=1e-6)


def test_synthetic_forward_parity_and_delta():
    sf = synthetic_forward(S, 100.0, T, R, SIG)
    assert sf.net.market_value == pytest.approx(S - math.exp(-R * T) * 100.0, abs=1e-9)
    assert sf.net.delta == pytest.approx(1.0, abs=1e-4)
    assert sf.net.gamma == pytest.approx(0.0, abs=1e-6)


def test_collar_matches_legs():
    c = collar(S, 90.0, 110.0, T, R, SIG)
    legs = price(S, 90.0, T, R, SIG, "put") - price(S, 110.0, T, R, SIG, "call")
    assert c.net.market_value == pytest.approx(legs, abs=1e-9)


def test_validation():
    with pytest.raises(ValueError):
        box_spread(S, 110.0, 90.0, T, R, SIG)
    with pytest.raises(ValueError):
        collar(S, 110.0, 90.0, T, R, SIG)
