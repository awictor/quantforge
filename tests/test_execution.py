"""Almgren-Chriss optimal execution."""

import pytest

import math

from quantforge import (
    execution_trajectory, execution_trades, expected_cost, cost_variance,
    efficient_frontier_point,
    kyle_lambda, kyle_impact, square_root_impact, implementation_shortfall,
    twap_schedule, vwap_schedule, pov_schedule,
)


X, N, T, SIG, ETA, GAMMA = 1e6, 10, 1.0, 0.3, 1e-6, 5e-7


def test_trajectory_endpoints():
    t = execution_trajectory(X, N, T, 0.0, SIG, ETA)
    assert t[0] == pytest.approx(X)
    assert t[-1] == pytest.approx(0.0, abs=1e-6)


def test_zero_lambda_is_linear():
    t = execution_trajectory(X, N, T, 0.0, SIG, ETA)
    assert all(t[k] == pytest.approx(X * (1 - k / N), abs=1e-3) for k in range(N + 1))


def test_trajectory_monotone_decreasing():
    t = execution_trajectory(X, N, T, 1e-6, SIG, ETA)
    assert all(t[k] >= t[k + 1] - 1e-6 for k in range(N))


def test_trades_sum_to_total():
    t = execution_trajectory(X, N, T, 0.0, SIG, ETA)
    assert sum(execution_trades(t)) == pytest.approx(X, abs=1e-3)


def test_risk_aversion_front_loads():
    linear = execution_trajectory(X, N, T, 0.0, SIG, ETA)
    urgent = execution_trajectory(X, N, T, 1e-6, SIG, ETA)
    assert all(urgent[k] <= linear[k] + 1e-3 for k in range(N + 1))
    assert urgent[1] < linear[1]


def test_efficient_frontier_tradeoff():
    c0, v0 = efficient_frontier_point(X, N, T, 0.0, SIG, ETA, GAMMA)
    c1, v1 = efficient_frontier_point(X, N, T, 1e-6, SIG, ETA, GAMMA)
    assert v1 < v0   # more urgency -> less timing risk
    assert c1 > c0   # ... at higher impact cost


def test_permanent_cost_component():
    # Two-point trajectory: all permanent + one temporary term.
    perm = 0.5 * GAMMA * X * X
    assert expected_cost([X, 0.0], T, GAMMA, 0.0) == pytest.approx(perm)


def test_variance_falls_with_front_loading():
    linear = execution_trajectory(X, N, T, 0.0, SIG, ETA)
    urgent = execution_trajectory(X, N, T, 1e-6, SIG, ETA)
    assert cost_variance(urgent, T, SIG) < cost_variance(linear, T, SIG)


def test_twap_equal_slices():
    t = twap_schedule(1e6, 10)
    assert all(x == pytest.approx(1e5) for x in t)
    assert sum(t) == pytest.approx(1e6)


def test_vwap_proportional_to_volume():
    v = vwap_schedule(1e6, [100, 200, 300, 400])
    assert v[1] / v[0] == pytest.approx(2.0)
    assert sum(v) == pytest.approx(1e6)


def test_vwap_flat_profile_is_twap():
    v = vwap_schedule(1e6, [1, 1, 1, 1])
    assert all(x == pytest.approx(2.5e5) for x in v)


def test_pov_proportional_to_interval_volume():
    prof = [100, 200, 300, 400]
    p = pov_schedule(prof, 0.1)
    assert all(p[i] == pytest.approx(0.1 * prof[i]) for i in range(4))


def test_pov_caps_at_total():
    p = pov_schedule([100, 200, 300, 400], 0.5, total_shares=120)
    assert sum(p) == pytest.approx(120)


def test_schedule_validation():
    with pytest.raises(ValueError):
        pov_schedule([100, 200], 1.5)
    with pytest.raises(ValueError):
        vwap_schedule(1e6, [0, 0, 0])
    with pytest.raises(ValueError):
        twap_schedule(1e6, 0)


def test_kyle_impact_linear():
    assert kyle_impact(2e5, 0.3, 1e7) == pytest.approx(2 * kyle_impact(1e5, 0.3, 1e7))


def test_kyle_lambda_formula_and_monotonicity():
    assert kyle_lambda(0.3, 1e7) == pytest.approx(0.3 / 1e7)
    assert kyle_lambda(0.6, 1e7) > kyle_lambda(0.3, 1e7)
    assert kyle_lambda(0.3, 2e7) < kyle_lambda(0.3, 1e7)


def test_square_root_law():
    # Quadrupling size doubles impact.
    assert square_root_impact(4e5, 0.3, 1e7) == pytest.approx(
        2 * square_root_impact(1e5, 0.3, 1e7))


def test_square_root_is_concave_per_share():
    assert square_root_impact(2e5, 0.3, 1e7) / 2e5 < \
        square_root_impact(1e5, 0.3, 1e7) / 1e5


def test_implementation_shortfall_decomposition():
    tr = execution_trajectory(X, N, T, 1e-6, SIG, ETA)
    perm, temp, tstd, tot = implementation_shortfall(tr, T, GAMMA, ETA, SIG)
    assert perm + temp == pytest.approx(tot, abs=1e-3)
    assert perm == pytest.approx(0.5 * GAMMA * X * X, abs=1e-3)
    assert tstd == pytest.approx(math.sqrt(cost_variance(tr, T, SIG)))


def test_impact_validation():
    with pytest.raises(ValueError):
        square_root_impact(1e5, 0.3, 0)
    with pytest.raises(ValueError):
        kyle_lambda(0.3, 0)


def test_validation():
    with pytest.raises(ValueError):
        execution_trajectory(-1, N, T, 0.0, SIG, ETA)
    with pytest.raises(ValueError):
        execution_trajectory(X, 0, T, 0.0, SIG, ETA)
    with pytest.raises(ValueError):
        execution_trajectory(X, N, T, 0.0, SIG, 0.0)
