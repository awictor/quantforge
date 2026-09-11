"""Commodity cost-of-carry forward pricing."""

import math

import pytest

from quantforge import (
    commodity_forward, implied_convenience_yield, implied_storage_cost,
    net_cost_of_carry, commodity_forward_curve, is_backwardation,
    commodity_calendar_spread, convenience_yield_curve, seasonal_forward,
    schwartz_log_mean, schwartz_log_variance, schwartz_forward,
    schwartz_option, mean_reversion_half_life, schwartz_implied_alpha,
    roll_yield, carry_roll_yield, schwartz_futures_volatility,
    margrabe_exchange_option, kirk_spread_option,
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


def test_schwartz_option_put_call_parity():
    F = schwartz_forward(SK, SKAPPA, SALPHA, SSIG, 2.0)
    c = schwartz_option(SK, SKAPPA, SALPHA, SSIG, 55.0, 0.05, 2.0, True)
    p = schwartz_option(SK, SKAPPA, SALPHA, SSIG, 55.0, 0.05, 2.0, False)
    assert c - p == pytest.approx(math.exp(-0.05 * 2.0) * (F - 55.0), abs=1e-9)


def test_schwartz_option_zero_vol_intrinsic():
    F = schwartz_forward(SK, SKAPPA, SALPHA, 0.0, 2.0)
    c = schwartz_option(SK, SKAPPA, SALPHA, 0.0, 55.0, 0.05, 2.0, True)
    assert c == pytest.approx(math.exp(-0.05 * 2.0) * max(F - 55.0, 0.0), abs=1e-9)


def test_mean_reversion_half_life():
    assert mean_reversion_half_life(SKAPPA) == pytest.approx(math.log(2) / SKAPPA)
    assert mean_reversion_half_life(3.0) < mean_reversion_half_life(1.0)


def test_schwartz_implied_alpha_round_trip():
    F = schwartz_forward(SK, SKAPPA, SALPHA, SSIG, 2.0)
    alpha = schwartz_implied_alpha(SK, F, SKAPPA, SSIG, 2.0)
    assert alpha == pytest.approx(SALPHA, abs=1e-9)
    assert schwartz_forward(SK, SKAPPA, alpha, SSIG, 2.0) == pytest.approx(F, abs=1e-9)


def test_schwartz_option_calibration_validation():
    with pytest.raises(ValueError):
        schwartz_implied_alpha(SK, 50.0, SKAPPA, SSIG, 0)
    with pytest.raises(ValueError):
        schwartz_option(SK, SKAPPA, SALPHA, SSIG, -1, 0.05, 2.0)
    with pytest.raises(ValueError):
        mean_reversion_half_life(0)


def test_roll_yield_positive_in_backwardation():
    Fn = commodity_forward(S, R, 1, 0.0, 0.10)
    Ff = commodity_forward(S, R, 2, 0.0, 0.10)
    assert roll_yield(Fn, Ff, 1, 2) > 0


def test_roll_yield_negative_in_contango():
    Fn = commodity_forward(S, R, 1, 0.02, 0.01)
    Ff = commodity_forward(S, R, 2, 0.02, 0.01)
    assert roll_yield(Fn, Ff, 1, 2) < 0


def test_carry_roll_yield_is_negative_net_carry():
    assert carry_roll_yield(R, 0.02, 0.03) == pytest.approx(
        -net_cost_of_carry(R, 0.02, 0.03))


def test_roll_yield_matches_carry_model():
    Fn = commodity_forward(S, R, 1, 0.02, 0.01)
    Ff = commodity_forward(S, R, 2, 0.02, 0.01)
    assert roll_yield(Fn, Ff, 1, 2) == pytest.approx(carry_roll_yield(R, 0.02, 0.01))


def test_schwartz_futures_vol_samuelson():
    assert schwartz_futures_volatility(0.3, 1.5, 0) == pytest.approx(0.3)
    vs = [schwartz_futures_volatility(0.3, 1.5, T) for T in (0, 0.5, 1, 2, 5)]
    assert all(vs[i] > vs[i + 1] for i in range(len(vs) - 1))


def test_roll_yield_validation():
    with pytest.raises(ValueError):
        roll_yield(-1, 100, 1, 2)
    with pytest.raises(ValueError):
        roll_yield(100, 100, 2, 1)
    with pytest.raises(ValueError):
        schwartz_futures_volatility(0.3, 0, 1)


SP_F1, SP_F2, SP_S1, SP_S2, SP_RHO, SP_R, SP_T = 100.0, 90.0, 0.3, 0.25, 0.4, 0.05, 1.0


def test_kirk_reduces_to_margrabe_at_zero_strike():
    m = margrabe_exchange_option(SP_F1, SP_F2, SP_S1, SP_S2, SP_RHO, SP_R, SP_T)
    k = kirk_spread_option(SP_F1, SP_F2, 0.0, SP_S1, SP_S2, SP_RHO, SP_R, SP_T, True)
    assert k == pytest.approx(m, abs=1e-9)


def test_kirk_put_call_parity():
    c = kirk_spread_option(SP_F1, SP_F2, 5.0, SP_S1, SP_S2, SP_RHO, SP_R, SP_T, True)
    p = kirk_spread_option(SP_F1, SP_F2, 5.0, SP_S1, SP_S2, SP_RHO, SP_R, SP_T, False)
    assert c - p == pytest.approx(math.exp(-SP_R * SP_T) * (SP_F1 - SP_F2 - 5.0), abs=1e-9)


def test_spread_option_correlation_monotone():
    lo = kirk_spread_option(SP_F1, SP_F2, 5.0, SP_S1, SP_S2, 0.0, SP_R, SP_T)
    hi = kirk_spread_option(SP_F1, SP_F2, 5.0, SP_S1, SP_S2, 0.8, SP_R, SP_T)
    assert hi < lo  # higher correlation shrinks the spread vol
    assert margrabe_exchange_option(SP_F1, SP_F2, SP_S1, SP_S2, 0.9, SP_R, SP_T) < \
        margrabe_exchange_option(SP_F1, SP_F2, SP_S1, SP_S2, 0.0, SP_R, SP_T)


def test_spread_option_positive():
    assert margrabe_exchange_option(SP_F1, SP_F2, SP_S1, SP_S2, SP_RHO, SP_R, SP_T) > 0
    assert kirk_spread_option(SP_F1, SP_F2, 5.0, SP_S1, SP_S2, SP_RHO, SP_R, SP_T) > 0


def test_spread_option_validation():
    with pytest.raises(ValueError):
        kirk_spread_option(SP_F1, SP_F2, -200, SP_S1, SP_S2, SP_RHO, SP_R, SP_T)
    with pytest.raises(ValueError):
        margrabe_exchange_option(-1, SP_F2, SP_S1, SP_S2, SP_RHO, SP_R, SP_T)


def test_validation():
    with pytest.raises(ValueError):
        implied_convenience_yield(S, 100, R, 0)
    with pytest.raises(ValueError):
        commodity_forward(-1, R, T)
    with pytest.raises(ValueError):
        commodity_forward(S, R, -1)
