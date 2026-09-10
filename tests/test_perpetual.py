"""Tests for perpetual American options."""

import math

import pytest

from quantforge import (
    perpetual_american, perpetual_exercise_boundary, american_price, OptionType,
)


@pytest.mark.slow
def test_perpetual_put_matches_long_maturity_tree():
    S, K, r, sigma = 100, 100, 0.08, 0.3
    perp = perpetual_american(S, K, r, sigma, OptionType.PUT)
    tree = american_price(S, K, 50.0, r, sigma, OptionType.PUT, b=r, steps=3000)
    assert perp == pytest.approx(tree, abs=0.1)


@pytest.mark.slow
def test_perpetual_call_with_dividend_matches_long_tree():
    S, K, r, sigma, b = 100, 100, 0.08, 0.3, 0.03
    perp = perpetual_american(S, K, r, sigma, OptionType.CALL, b=b)
    tree = american_price(S, K, 50.0, r, sigma, OptionType.CALL, b=b, steps=3000)
    assert perp == pytest.approx(tree, abs=0.2)


def test_put_dominates_intrinsic():
    v = perpetual_american(80, 100, 0.08, 0.3, OptionType.PUT)
    assert v >= (100 - 80) - 1e-9


def test_call_dominates_intrinsic():
    v = perpetual_american(150, 100, 0.08, 0.3, OptionType.CALL, b=0.03)
    assert v >= (150 - 100) - 1e-9


def test_deep_itm_put_near_intrinsic_beyond_boundary():
    # Below the exercise boundary the put equals its intrinsic value.
    K, r, sigma = 100, 0.08, 0.3
    Sb = perpetual_exercise_boundary(K, r, sigma, OptionType.PUT)
    S = 0.5 * Sb  # well past the boundary
    assert perpetual_american(S, K, r, sigma, OptionType.PUT) == pytest.approx(K - S)


def test_call_never_exercised_when_carry_ge_rate():
    # b >= r: a perpetual call is never exercised; its value tends to the spot.
    v = perpetual_american(100, 100, 0.05, 0.3, OptionType.CALL, b=0.05)
    assert v == pytest.approx(100.0)


def test_boundary_put_below_strike_call_above():
    K, r, sigma = 100, 0.08, 0.3
    put_b = perpetual_exercise_boundary(K, r, sigma, OptionType.PUT)
    call_b = perpetual_exercise_boundary(K, r, sigma, OptionType.CALL, b=0.03)
    assert put_b < K < call_b


def test_higher_vol_raises_value():
    lo = perpetual_american(100, 100, 0.08, 0.2, OptionType.PUT)
    hi = perpetual_american(100, 100, 0.08, 0.4, OptionType.PUT)
    assert hi > lo


def test_rejects_bad_inputs():
    with pytest.raises(ValueError):
        perpetual_american(100, 100, 0.0, 0.3, OptionType.PUT)   # r must be > 0
    with pytest.raises(ValueError):
        perpetual_american(100, 100, 0.08, 0.0, OptionType.PUT)  # sigma > 0
