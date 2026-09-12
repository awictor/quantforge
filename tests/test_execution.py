"""Almgren-Chriss optimal execution."""

import pytest

from quantforge import (
    execution_trajectory, execution_trades, expected_cost, cost_variance,
    efficient_frontier_point,
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


def test_validation():
    with pytest.raises(ValueError):
        execution_trajectory(-1, N, T, 0.0, SIG, ETA)
    with pytest.raises(ValueError):
        execution_trajectory(X, 0, T, 0.0, SIG, ETA)
    with pytest.raises(ValueError):
        execution_trajectory(X, N, T, 0.0, SIG, 0.0)
