"""Risk-of-ruin and first-passage drawdown probabilities."""

import math
import random

import pytest

from quantforge import (gamblers_ruin_probability, risk_of_ruin_units,
                        ruin_probability_gbm)


def _recursion_ruin(start, target, p):
    """Solve P[i] = p P[i+1] + (1-p) P[i-1], P[0]=1, P[N]=0 by iteration."""
    n = target
    prob = [0.0] * (n + 1)
    prob[0] = 1.0
    for _ in range(200000):
        new = prob[:]
        for i in range(1, n):
            new[i] = p * prob[i + 1] + (1 - p) * prob[i - 1]
        if max(abs(new[i] - prob[i]) for i in range(n + 1)) < 1e-14:
            prob = new
            break
        prob = new
    return prob[start]


def test_gamblers_ruin_matches_recursion():
    for p in (0.45, 0.5, 0.55):
        assert abs(gamblers_ruin_probability(3, 10, p) - _recursion_ruin(3, 10, p)) < 1e-6


def test_fair_game_is_linear():
    assert abs(gamblers_ruin_probability(4, 10, 0.5) - 0.6) < 1e-12
    assert abs(gamblers_ruin_probability(1, 5, 0.5) - 0.8) < 1e-12


def test_favorable_game_lower_ruin():
    fav = gamblers_ruin_probability(5, 20, 0.55)
    unfav = gamblers_ruin_probability(5, 20, 0.45)
    assert fav < 0.5 < unfav


def test_risk_of_ruin_units_falls_with_bankroll():
    a = risk_of_ruin_units(0.55, 5)
    b = risk_of_ruin_units(0.55, 20)
    assert b < a


def test_ruin_gbm_matches_first_passage_formula():
    mu, sigma, loss = 0.15, 0.30, 0.30
    a = -math.log(1 - loss)
    closed = math.exp(-2 * mu * a / (sigma * sigma))
    assert abs(ruin_probability_gbm(loss, mu, sigma) - closed) < 1e-12


def test_ruin_gbm_monte_carlo():
    mu, sigma, loss = 0.20, 0.25, 0.25
    a = -math.log(1 - loss)
    rng = random.Random(31)
    nsim, steps, horizon = 12000, 5000, 300.0
    dt = horizon / steps
    sq = math.sqrt(dt)
    hit = 0
    for _ in range(nsim):
        x = 0.0
        for _ in range(steps):
            x += mu * dt + sigma * sq * rng.gauss(0, 1)
            if x <= -a:
                hit += 1
                break
    # Finite horizon truncates late hits, so MC is a lower bound on the closed form.
    mc = hit / nsim
    closed = ruin_probability_gbm(loss, mu, sigma)
    assert mc <= closed + 0.02
    assert closed - mc < 0.05


def test_ruin_gbm_decreasing_in_edge_and_one_without_drift():
    weak = ruin_probability_gbm(0.3, 0.05, 0.3)
    strong = ruin_probability_gbm(0.3, 0.20, 0.3)
    assert strong < weak
    assert ruin_probability_gbm(0.3, 0.0, 0.3) == 1.0
    assert ruin_probability_gbm(0.3, -0.1, 0.3) == 1.0


def test_validation():
    with pytest.raises(ValueError):
        gamblers_ruin_probability(10, 10, 0.5)     # start not < target
    with pytest.raises(ValueError):
        gamblers_ruin_probability(3, 10, 1.0)
    with pytest.raises(ValueError):
        risk_of_ruin_units(0.55, 0)
    with pytest.raises(ValueError):
        ruin_probability_gbm(1.0, 0.1, 0.2)
    with pytest.raises(ValueError):
        ruin_probability_gbm(0.3, 0.1, 0.0)
