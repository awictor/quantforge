"""Tests for the double-knockout (corridor) barrier Monte Carlo pricer."""

import pytest

from quantforge import double_knockout_mc, call_price, OptionType


@pytest.mark.slow
def test_wide_barriers_approach_vanilla():
    S, K = 100, 100
    v = double_knockout_mc(S, K, 1.0, 0.05, 0.25, lower=1.0, upper=1e6,
                           option_type=OptionType.CALL, n_steps=300,
                           n_paths=60_000, seed=1)
    assert v.price == pytest.approx(call_price(S, K, 1.0, 0.05, 0.25), abs=0.1)


@pytest.mark.slow
def test_tighter_corridor_is_cheaper():
    wide = double_knockout_mc(100, 100, 1.0, 0.05, 0.25, lower=70, upper=140,
                              n_steps=200, n_paths=40_000, seed=2)
    tight = double_knockout_mc(100, 100, 1.0, 0.05, 0.25, lower=90, upper=115,
                               n_steps=200, n_paths=40_000, seed=2)
    assert tight.price < wide.price


@pytest.mark.slow
def test_rebate_raises_value():
    no_reb = double_knockout_mc(100, 100, 1.0, 0.05, 0.3, lower=90, upper=115,
                                rebate=0.0, n_steps=200, n_paths=40_000, seed=3)
    with_reb = double_knockout_mc(100, 100, 1.0, 0.05, 0.3, lower=90, upper=115,
                                  rebate=5.0, n_steps=200, n_paths=40_000, seed=3)
    assert with_reb.price > no_reb.price


def test_price_nonnegative():
    v = double_knockout_mc(100, 100, 1.0, 0.05, 0.25, lower=85, upper=120,
                           n_steps=100, n_paths=20_000, seed=4)
    assert v.price >= 0


def test_reproducible_with_seed():
    a = double_knockout_mc(100, 100, 0.5, 0.05, 0.25, lower=90, upper=115,
                           n_steps=50, n_paths=5_000, seed=9)
    b = double_knockout_mc(100, 100, 0.5, 0.05, 0.25, lower=90, upper=115,
                           n_steps=50, n_paths=5_000, seed=9)
    assert a.price == b.price


def test_rejects_spot_outside_corridor():
    with pytest.raises(ValueError):
        double_knockout_mc(100, 100, 1.0, 0.05, 0.25, lower=105, upper=120, n_paths=100)
