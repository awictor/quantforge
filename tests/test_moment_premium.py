"""Skewness- and kurtosis-risk premia: realized vs BKM-implied."""

import math
import random

import pytest

from quantforge import moment_risk_premia


def _gbm(sig=0.2, n=252, seed=1):
    rng = random.Random(seed)
    closes = [100.0]
    for _ in range(n):
        closes.append(closes[-1] * math.exp(-0.5 * sig * sig / 252
                                            + sig / math.sqrt(252) * rng.gauss(0, 1)))
    return closes


def _down(K):
    return max(0.05, 0.2 + 0.2 * math.log(100.0 / K))


def _flat(K):
    return 0.2


def test_downward_smile_negative_skew_premium():
    res = moment_risk_premia(_gbm(), 100, 0.5, 0.03, _down)
    # A downward-skewed smile implies more-negative skew than a near-symmetric
    # realized series -> negative skew premium (crash-protection demand).
    assert res["implied_skew"] < res["realized_skew"]
    assert res["skew_premium"] < 0.0


def test_downward_smile_positive_kurt_premium():
    res = moment_risk_premia(_gbm(), 100, 0.5, 0.03, _down)
    assert res["kurt_premium"] > 0.0


def test_flat_smile_small_skew_premium():
    res = moment_risk_premia(_gbm(), 100, 0.5, 0.03, _flat)
    # Both realized (GBM) and implied (flat) skew are near zero.
    assert abs(res["skew_premium"]) < 0.2


def test_premium_components_consistent():
    res = moment_risk_premia(_gbm(), 100, 0.5, 0.03, _down)
    assert res["skew_premium"] == pytest.approx(
        res["implied_skew"] - res["realized_skew"])
    assert res["kurt_premium"] == pytest.approx(
        res["implied_kurt"] - res["realized_kurt"])


def test_too_few_closes_raises():
    with pytest.raises(ValueError):
        moment_risk_premia([100, 101], 100, 0.5, 0.03, _flat)
