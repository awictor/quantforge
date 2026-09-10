"""Tests for vega term-structure bucketing."""

import pytest

from quantforge import Contract, vega_buckets, price_book, vega, OptionType


def _leg(t, qty=1, K=100):
    return Contract(S=100, K=K, t=t, r=0.04, sigma=0.25, option_type="call",
                    qty=qty, multiplier=1)


def test_buckets_sum_to_net_vega():
    book = [_leg(0.1), _leg(0.4), _leg(1.5), _leg(3.0), _leg(7.0)]
    vb = vega_buckets(book)
    net = price_book(book).net.vega
    assert vb.total == pytest.approx(net, abs=1e-9)


def test_contract_lands_in_correct_bucket():
    # A single 9-month option with default edges (…, 1.0y, …) -> "0.5-1y".
    vb = vega_buckets([_leg(0.75)], edges=(0.25, 0.5, 1.0, 2.0, 5.0))
    assert vb.buckets["0.5-1y"] == pytest.approx(vega(100, 100, 0.75, 0.04, 0.25))
    # Every other bucket is empty.
    for lab, val in vb.buckets.items():
        if lab != "0.5-1y":
            assert val == pytest.approx(0.0)


def test_long_expiry_goes_to_beyond_bucket():
    vb = vega_buckets([_leg(8.0)], edges=(0.25, 0.5, 1.0, 2.0, 5.0))
    assert vb.buckets[">5y"] == pytest.approx(vega(100, 100, 8.0, 0.04, 0.25))


def test_front_vs_back_separation():
    # Long front-month vol, short back-month.
    book = [_leg(0.2, qty=5), _leg(3.0, qty=-5)]
    vb = vega_buckets(book)
    assert vb.buckets["<=0.25y"] > 0
    assert vb.buckets["2-5y"] < 0


def test_short_position_negative_vega():
    vb = vega_buckets([_leg(1.5, qty=-3)])
    assert vb.buckets["1-2y"] < 0


def test_boundary_expiry_goes_to_lower_bucket():
    # An expiry exactly on an edge goes to the bucket whose edge equals it.
    vb = vega_buckets([_leg(1.0)], edges=(0.5, 1.0, 2.0))
    assert vb.buckets["0.5-1y"] == pytest.approx(vega(100, 100, 1.0, 0.04, 0.25))


def test_rejects_unsorted_edges():
    with pytest.raises(ValueError):
        vega_buckets([_leg(1.0)], edges=(1.0, 0.5))
