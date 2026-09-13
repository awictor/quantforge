"""Trade-sign classification: tick rule, quote rule, Lee-Ready."""

import pytest

from quantforge import tick_rule, quote_rule, lee_ready


def test_tick_rule_directions():
    # default +1, uptick, zero-tick carries, downtick, uptick
    assert tick_rule([10, 11, 11, 10, 12]) == [1, 1, 1, -1, 1]


def test_tick_rule_all_flat_carry_default():
    assert tick_rule([5, 5, 5]) == [1, 1, 1]


def test_quote_rule_above_below_at_mid():
    prices = [10.3, 9.7, 10.0]
    bids = [10.0, 9.5, 9.9]
    asks = [10.4, 10.1, 10.1]        # mids 10.2, 9.8, 10.0
    assert quote_rule(prices, bids, asks) == [1, -1, 0]


def test_lee_ready_resolves_midpoint_with_tick():
    prices = [10.3, 9.7, 10.0]
    bids = [10.0, 9.5, 9.9]
    asks = [10.4, 10.1, 10.1]
    # third trade at the mid -> tick rule: 10.0 > 9.7 -> +1
    assert lee_ready(prices, bids, asks) == [1, -1, 1]


def test_lee_ready_has_no_zeros():
    prices = [10.0, 10.0, 10.0]      # all at mid
    bids = [9.9, 9.9, 9.9]
    asks = [10.1, 10.1, 10.1]
    signs = lee_ready(prices, bids, asks)
    assert all(s in (-1, 1) for s in signs)


def test_feeds_order_flow_measures():
    from quantforge import vpin, order_flow_imbalance
    prices = [10.3, 9.7, 10.0, 10.5, 9.6]
    bids = [10, 9.5, 9.9, 10.2, 9.4]
    asks = [10.4, 10.1, 10.1, 10.6, 9.9]
    vols = [100, 200, 150, 300, 250]
    signs = lee_ready(prices, bids, asks)
    buys = [vols[i] if signs[i] > 0 else 0 for i in range(len(vols))]
    sells = [vols[i] if signs[i] < 0 else 0 for i in range(len(vols))]
    assert 0.0 <= vpin(buys, sells) <= 1.0
    assert -1.0 <= order_flow_imbalance(buys, sells) <= 1.0


def test_validation():
    with pytest.raises(ValueError):
        tick_rule([])
    with pytest.raises(ValueError):
        quote_rule([10.0], [9.9], [])
    with pytest.raises(ValueError):
        lee_ready([10.0, 10.1], [9.9], [10.1])
