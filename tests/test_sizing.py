"""Tests for position-sizing / hedge-quantity helpers."""

import pytest

from quantforge import (
    Contract, price_book, greeks,
    delta_hedge_shares, neutralize, vega_neutral_quantity, gamma_neutral_quantity,
    OptionType,
)


def _book():
    return price_book([
        Contract(S=100, K=100, t=0.5, r=0.04, sigma=0.25, option_type="call",
                 qty=10, multiplier=1),
    ])


def test_delta_hedge_zeros_net_delta():
    book = _book()
    shares = delta_hedge_shares(book)
    # A share has delta 1; net delta + shares should be zero.
    assert book.net.delta + shares == pytest.approx(0.0, abs=1e-9)
    assert shares < 0  # long calls -> sell shares


def test_vega_neutral_zeros_net_vega():
    book = _book()
    # Hedge with a different-strike option.
    hk = dict(S=100, K=110, t=0.5, r=0.04, sigma=0.25, option_type=OptionType.CALL)
    q = vega_neutral_quantity(book, **hk)
    g = greeks(100, 110, 0.5, 0.04, 0.25, OptionType.CALL)
    assert book.net.vega + q * g.vega == pytest.approx(0.0, abs=1e-6)


def test_gamma_neutral_zeros_net_gamma():
    book = _book()
    hk = dict(S=100, K=95, t=0.5, r=0.04, sigma=0.25, option_type=OptionType.PUT)
    q = gamma_neutral_quantity(book, **hk)
    g = greeks(100, 95, 0.5, 0.04, 0.25, OptionType.PUT)
    assert book.net.gamma + q * g.gamma == pytest.approx(0.0, abs=1e-8)


def test_neutralize_hits_arbitrary_target():
    book = _book()
    hk = dict(S=100, K=100, t=0.5, r=0.04, sigma=0.25, option_type=OptionType.CALL)
    target = 3.0
    q = neutralize(book, "delta", target=target, **hk)
    g = greeks(100, 100, 0.5, 0.04, 0.25, OptionType.CALL)
    assert book.net.delta + q * g.delta == pytest.approx(target, abs=1e-6)


def test_neutralize_respects_multiplier():
    book = _book()
    hk = dict(S=100, K=100, t=0.5, r=0.04, sigma=0.25, option_type=OptionType.CALL)
    q = neutralize(book, "vega", multiplier=100, **hk)
    g = greeks(100, 100, 0.5, 0.04, 0.25, OptionType.CALL)
    assert book.net.vega + q * 100 * g.vega == pytest.approx(0.0, abs=1e-6)


def test_neutralize_rejects_bad_greek():
    with pytest.raises(ValueError):
        neutralize(_book(), "rho", S=100, K=100, t=0.5, r=0.04, sigma=0.25)


def test_neutralize_rejects_zero_greek_hedge():
    # A hedge at t=0 has zero vega, so vega-neutralizing is impossible.
    with pytest.raises(ValueError):
        vega_neutral_quantity(_book(), S=100, K=100, t=0.0, r=0.04, sigma=0.25)
