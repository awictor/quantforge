"""Commodity cost-of-carry forward pricing."""

import math

import pytest

from quantforge import (
    commodity_forward, implied_convenience_yield, implied_storage_cost,
    net_cost_of_carry, commodity_forward_curve, is_backwardation,
)


S, R, T, U, Y = 100.0, 0.05, 2.0, 0.02, 0.03


def test_carry_parity():
    assert commodity_forward(S, R, T, U, Y) == pytest.approx(
        S * math.exp((R + U - Y) * T), abs=1e-9)


def test_implied_convenience_yield_inverts():
    F = commodity_forward(S, R, T, U, Y)
    assert implied_convenience_yield(S, F, R, T, U) == pytest.approx(Y, abs=1e-12)


def test_implied_storage_cost_inverts():
    F = commodity_forward(S, R, T, U, Y)
    assert implied_storage_cost(S, F, R, T, Y) == pytest.approx(U, abs=1e-12)


def test_net_cost_of_carry():
    assert net_cost_of_carry(R, U, Y) == pytest.approx(R + U - Y)


def test_contango_forward_above_spot():
    assert commodity_forward(S, 0.05, T, 0.02, 0.01) > S
    assert not is_backwardation(0.05, 0.02, 0.01)


def test_backwardation_forward_below_spot():
    assert commodity_forward(S, 0.05, T, 0.0, 0.10) < S
    assert is_backwardation(0.05, 0.0, 0.10)


def test_forward_curve_rises_in_contango():
    curve = commodity_forward_curve(S, 0.05, [0.5, 1, 2, 5], 0.02, 0.01)
    assert all(curve[i][1] < curve[i + 1][1] for i in range(len(curve) - 1))


def test_forward_at_zero_maturity_is_spot():
    assert commodity_forward(S, R, 0.0, U, Y) == pytest.approx(S)


def test_validation():
    with pytest.raises(ValueError):
        implied_convenience_yield(S, 100, R, 0)
    with pytest.raises(ValueError):
        commodity_forward(-1, R, T)
    with pytest.raises(ValueError):
        commodity_forward(S, R, -1)
