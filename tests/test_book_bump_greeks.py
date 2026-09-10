"""Tests for whole-book bump-and-reprice Greeks."""

import pytest

from quantforge import Contract, book_bump_greeks, price_book, greeks, OptionType


def _book():
    return [
        Contract(S=100, K=100, t=0.5, r=0.04, sigma=0.25, option_type="call", qty=10),
        Contract(S=100, K=95, t=0.5, r=0.04, sigma=0.25, option_type="put", qty=-5),
    ]


def test_matches_analytic_net_greeks():
    b = book_bump_greeks(_book())
    net = price_book(_book()).net
    assert b.delta == pytest.approx(net.delta, abs=1e-3)
    assert b.gamma == pytest.approx(net.gamma, abs=1e-4)
    assert b.vega == pytest.approx(net.vega, abs=1e-1)
    assert b.theta == pytest.approx(net.theta, abs=1e-1)


def test_price_matches_book_value():
    b = book_bump_greeks(_book())
    assert b.price == pytest.approx(price_book(_book()).net.market_value, abs=1e-9)


def test_single_leg_matches_scalar_delta():
    leg = [Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2, option_type="call", qty=1)]
    b = book_bump_greeks(leg)
    assert b.delta == pytest.approx(greeks(100, 100, 1.0, 0.05, 0.2, OptionType.CALL).delta,
                                    abs=1e-4)


def test_short_book_flips_delta():
    long = book_bump_greeks([Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                                      option_type="call", qty=1)])
    short = book_bump_greeks([Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                                       option_type="call", qty=-1)])
    assert short.delta == pytest.approx(-long.delta, abs=1e-6)


def test_empty_book_rejected():
    with pytest.raises(ValueError):
        book_bump_greeks([])
