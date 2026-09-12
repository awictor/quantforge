"""Augmented Dickey-Fuller and Engle-Granger cointegration tests."""

import random

import pytest

from quantforge import adf_test, engle_granger


def _ar1(n, phi, seed, sigma=1.0):
    rng = random.Random(seed)
    y = [0.0]
    for _ in range(n):
        y.append(phi * y[-1] + rng.gauss(0, sigma))
    return y[1:]


def _walk(n, seed):
    rng = random.Random(seed)
    y = [0.0]
    for _ in range(n):
        y.append(y[-1] + rng.gauss(0, 1))
    return y


def test_stationary_series_rejects_unit_root():
    r = adf_test(_ar1(1000, 0.5, 1))
    assert r["reject_5pct"]
    assert r["statistic"] < -2.86


def test_random_walk_fails_to_reject():
    # A single ADF draw has a ~5% false-positive rate by construction, so require
    # the test to hold on the clear majority of independent random walks.
    no_reject = sum(not adf_test(_walk(1000, s))["reject_5pct"] for s in range(20))
    assert no_reject >= 17


def test_cointegrated_pair_detected():
    x = _walk(1000, 3)
    noise = _ar1(len(x) + 1, 0.5, 4, sigma=0.3)[: len(x)]
    y = [2.0 * x[i] + noise[i] for i in range(len(x))]
    c = engle_granger(y, x)
    assert abs(c["hedge_ratio"] - 2.0) < 0.1
    assert c["cointegrated_5pct"]


def test_independent_walks_not_cointegrated():
    x = _walk(1000, 5)
    y = _walk(1000, 6)
    c = engle_granger(y, x)
    assert not c["cointegrated_5pct"]


def test_augmented_lags_run():
    r = adf_test(_ar1(1000, 0.5, 7), lags=4)
    assert r["reject_5pct"]


def test_validation():
    with pytest.raises(ValueError):
        adf_test([1.0, 2.0, 3.0])                 # too short
    with pytest.raises(ValueError):
        adf_test([1.0] * 20, lags=-1)             # negative lags
    with pytest.raises(ValueError):
        engle_granger([1.0, 2.0], [1.0])          # length mismatch
