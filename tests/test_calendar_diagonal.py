"""Calendar and diagonal spread builders (strategy module)."""

import pytest

from quantforge import calendar_spread, diagonal_spread
from quantforge.bsm import price


S, K, R, SIG = 100.0, 100.0, 0.05, 0.2


def test_calendar_net_is_far_minus_near():
    b = calendar_spread(S, K, 0.25, 1.0, R, SIG, "call")
    near = price(S, K, 0.25, R, SIG, "call")
    far = price(S, K, 1.0, R, SIG, "call")
    assert b.net.market_value == pytest.approx(far - near, abs=1e-9)


def test_long_calendar_is_debit_and_long_vega():
    b = calendar_spread(S, K, 0.25, 1.0, R, SIG, "call")
    assert b.net.market_value > 0.0     # net debit (far worth more)
    assert b.net.vega > 0.0             # long the longer-dated vega


def test_diagonal_net_matches_legs():
    d = diagonal_spread(S, 95.0, 105.0, 0.25, 1.0, R, SIG, "call")
    near = price(S, 95.0, 0.25, R, SIG, "call")
    far = price(S, 105.0, 1.0, R, SIG, "call")
    assert d.net.market_value == pytest.approx(far - near, abs=1e-9)


def test_put_calendar():
    b = calendar_spread(S, K, 0.25, 1.0, R, SIG, "put")
    near = price(S, K, 0.25, R, SIG, "put")
    far = price(S, K, 1.0, R, SIG, "put")
    assert b.net.market_value == pytest.approx(far - near, abs=1e-9)


def test_expiry_order_validation():
    with pytest.raises(ValueError):
        calendar_spread(S, K, 1.0, 0.25, R, SIG, "call")
    with pytest.raises(ValueError):
        diagonal_spread(S, 95.0, 105.0, 1.0, 0.25, R, SIG, "call")
