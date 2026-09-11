"""Commodity cost-of-carry forward pricing."""

import math

import pytest

from quantforge import (
    commodity_forward, implied_convenience_yield, implied_storage_cost,
    net_cost_of_carry, commodity_forward_curve, is_backwardation,
    commodity_calendar_spread, convenience_yield_curve, seasonal_forward,
    schwartz_log_mean, schwartz_log_variance, schwartz_forward,
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


def test_calendar_spread_sign_tracks_carry():
    assert commodity_calendar_spread(S, R, 1, 2, 0.02, 0.01) > 0   # contango
    assert commodity_calendar_spread(S, R, 1, 2, 0.0, 0.10) < 0    # backwardation


def test_convenience_yield_curve_reprices_forwards():
    quotes = commodity_forward_curve(S, R, [0.5, 1, 2], 0.02, 0.03)
    ys = convenience_yield_curve(S, R, quotes, 0.02)
    for i, (T, F) in enumerate(quotes):
        assert commodity_forward(S, R, T, 0.02, ys[i][1]) == pytest.approx(F, abs=1e-9)
        assert ys[i][1] == pytest.approx(0.03, abs=1e-12)


def test_seasonal_forward_scales():
    base = commodity_forward(S, R, 1, 0.02, 0.01)
    assert seasonal_forward(S, R, 1, 1.0, 0.02, 0.01) == pytest.approx(base)
    assert seasonal_forward(S, R, 1, 1.1, 0.02, 0.01) > base


def test_spread_seasonal_validation():
    with pytest.raises(ValueError):
        commodity_calendar_spread(S, R, 2, 1)
    with pytest.raises(ValueError):
        seasonal_forward(S, R, 1, -1)


SK, SKAPPA, SALPHA, SSIG = 50.0, 1.5, 4.0, 0.3


def test_schwartz_zero_maturity_is_spot():
    assert schwartz_log_mean(SK, SKAPPA, SALPHA, 0) == pytest.approx(math.log(SK))
    assert schwartz_log_variance(SSIG, SKAPPA, 0) == pytest.approx(0.0, abs=1e-12)
    assert schwartz_forward(SK, SKAPPA, SALPHA, SSIG, 0) == pytest.approx(SK, abs=1e-9)


def test_schwartz_long_run_limits():
    assert schwartz_log_mean(SK, SKAPPA, SALPHA, 100) == pytest.approx(SALPHA, abs=1e-9)
    assert schwartz_log_variance(SSIG, SKAPPA, 100) == pytest.approx(
        SSIG ** 2 / (2 * SKAPPA), abs=1e-9)
    assert schwartz_forward(SK, SKAPPA, SALPHA, SSIG, 100) == pytest.approx(
        math.exp(SALPHA + SSIG ** 2 / (4 * SKAPPA)), abs=1e-6)


def test_schwartz_variance_monotone():
    vs = [schwartz_log_variance(SSIG, SKAPPA, T) for T in (0.1, 0.5, 1, 2, 5, 10)]
    assert all(vs[i] < vs[i + 1] for i in range(len(vs) - 1))


def test_schwartz_mean_reverts_upward():
    # Spot below long-run level -> mean rises toward alpha.
    ms = [schwartz_log_mean(SK, SKAPPA, SALPHA, T) for T in (0.1, 0.5, 1, 2, 5)]
    assert all(ms[i] < ms[i + 1] for i in range(len(ms) - 1))


def test_schwartz_validation():
    with pytest.raises(ValueError):
        schwartz_forward(-1, SKAPPA, SALPHA, SSIG, 1)
    with pytest.raises(ValueError):
        schwartz_log_variance(SSIG, 0, 1)


def test_validation():
    with pytest.raises(ValueError):
        implied_convenience_yield(S, 100, R, 0)
    with pytest.raises(ValueError):
        commodity_forward(-1, R, T)
    with pytest.raises(ValueError):
        commodity_forward(S, R, -1)
