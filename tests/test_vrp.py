"""Variance risk premium: realized vs implied variance."""

import math
import random

import pytest

from quantforge import realized_variance, variance_risk_premium


def _gbm(sig, n=252, seed=1):
    rng = random.Random(seed)
    dt = 1.0 / 252.0
    closes = [100.0]
    for _ in range(n):
        closes.append(closes[-1] * math.exp((-0.5 * sig * sig) * dt
                                             + sig * math.sqrt(dt) * rng.gauss(0, 1)))
    return closes


def test_realized_variance_recovers_sigma_squared():
    rv = realized_variance(_gbm(0.2))
    assert rv == pytest.approx(0.04, abs=1e-2)


def test_negative_vrp_when_implied_above_realized():
    closes = _gbm(0.2)
    res = variance_risk_premium(closes, 0.06)
    assert res["vrp"] < 0.0
    assert res["vol_premium"] > 0.0
    assert res["ratio"] < 1.0


def test_positive_vrp_when_realized_above_implied():
    closes = _gbm(0.3)                       # realized ~0.09
    res = variance_risk_premium(closes, 0.04)
    assert res["vrp"] > 0.0
    assert res["ratio"] > 1.0


def test_components_consistent():
    closes = _gbm(0.25)
    res = variance_risk_premium(closes, 0.05)
    assert res["vrp"] == pytest.approx(res["realized_variance"] - 0.05)
    assert res["vol_premium"] == pytest.approx(
        math.sqrt(0.05) - math.sqrt(res["realized_variance"]))


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        realized_variance([100.0])
    with pytest.raises(ValueError):
        variance_risk_premium(_gbm(0.2), -0.01)
