"""Tests for book-level aggregate second-order Greeks."""

import pytest

from quantforge import (
    Contract, book_second_order, vanna, vomma, charm, veta, speed, zomma, color,
    OptionType,
)


def test_matches_scaled_leg_sum():
    contracts = [
        Contract(S=100, K=105, t=0.5, r=0.04, sigma=0.25, option_type="call",
                 qty=10, multiplier=100),
        Contract(S=100, K=95, t=0.5, r=0.04, sigma=0.30, option_type="put",
                 qty=-5, multiplier=100),
    ]
    bk = book_second_order(contracts)

    exp_vanna = (10 * 100 * vanna(100, 105, 0.5, 0.04, 0.25)
                 + -5 * 100 * vanna(100, 95, 0.5, 0.04, 0.30))
    assert bk.vanna == pytest.approx(exp_vanna, abs=1e-9)

    exp_charm = (10 * 100 * charm(100, 105, 0.5, 0.04, 0.25, OptionType.CALL)
                 + -5 * 100 * charm(100, 95, 0.5, 0.04, 0.30, OptionType.PUT))
    assert bk.charm == pytest.approx(exp_charm, abs=1e-9)


def test_all_fields_populated():
    bk = book_second_order([Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                                     option_type="call", qty=1)])
    for f in ("vanna", "vomma", "charm", "veta", "speed", "zomma", "color"):
        assert hasattr(bk, f)


def test_single_leg_equals_scalar():
    c = Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2, option_type="call", qty=1)
    bk = book_second_order([c])
    assert bk.vomma == pytest.approx(vomma(100, 100, 1.0, 0.05, 0.2))
    assert bk.speed == pytest.approx(speed(100, 100, 1.0, 0.05, 0.2))
    assert bk.zomma == pytest.approx(zomma(100, 100, 1.0, 0.05, 0.2))
    assert bk.color == pytest.approx(color(100, 100, 1.0, 0.05, 0.2))
    assert bk.veta == pytest.approx(veta(100, 100, 1.0, 0.05, 0.2))


def test_short_position_flips_sign():
    long = book_second_order([Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                                       option_type="call", qty=1)])
    short = book_second_order([Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                                        option_type="call", qty=-1)])
    assert short.vanna == pytest.approx(-long.vanna)
    assert short.vomma == pytest.approx(-long.vomma)


def test_vomma_is_nonnegative_for_long_options():
    # Long options have non-negative vomma (vol convexity) at typical strikes.
    for K in (80, 100, 120):
        bk = book_second_order([Contract(S=100, K=K, t=1.0, r=0.05, sigma=0.2,
                                         option_type="call", qty=1)])
        assert bk.vomma >= -1e-9


def test_offsetting_book_nets_to_zero():
    # Long and short the identical contract nets every second-order Greek to 0.
    legs = [
        Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2, option_type="call", qty=3),
        Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2, option_type="call", qty=-3),
    ]
    bk = book_second_order(legs)
    for f in ("vanna", "vomma", "charm", "veta", "speed", "zomma", "color"):
        assert getattr(bk, f) == pytest.approx(0.0, abs=1e-9)
