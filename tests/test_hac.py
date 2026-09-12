"""Newey-West HAC long-run variance."""

import math
import random

import pytest

from quantforge import (
    autocovariance, autocorrelation, newey_west_variance, newey_west_mean_se,
)


def _ar1(n, phi, seed, sigma=1.0):
    rng = random.Random(seed)
    y = [0.0]
    for _ in range(n + 1000):
        y.append(phi * y[-1] + rng.gauss(0, sigma))
    return y[1001:]


def test_lag_zero_is_sample_variance():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(500)]
    m = sum(x) / len(x)
    var = sum((v - m) ** 2 for v in x) / len(x)
    assert abs(newey_west_variance(x, 0) - var) < 1e-12


def test_always_nonnegative():
    rng = random.Random(2)
    x = [rng.gauss(0, 1) for _ in range(300)]
    assert all(newey_west_variance(x, L) >= 0.0 for L in (0, 1, 5, 20, 50))


def test_white_noise_near_sample_variance():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(500)]
    var = autocovariance(x, 0)
    assert abs(newey_west_variance(x, 10) - var) / var < 0.2


def test_ar1_approaches_long_run_variance():
    phi = 0.6
    y = _ar1(200000, phi, 5)
    lrv_true = 1.0 / (1 - phi) ** 2          # sigma_eps^2 / (1-phi)^2
    # Newey-West with a growing lag climbs toward the true long-run variance.
    nw5 = newey_west_variance(y, 5)
    nw50 = newey_west_variance(y, 50)
    assert nw5 < nw50 <= lrv_true * 1.05
    assert abs(nw50 - lrv_true) / lrv_true < 0.2


def test_autocorrelation_recovers_phi():
    phi = 0.6
    y = _ar1(200000, phi, 7)
    assert abs(autocorrelation(y, 1) - phi) < 0.02


def test_hac_se_exceeds_iid_for_positive_autocorrelation():
    phi = 0.6
    y = _ar1(50000, phi, 9)
    iid_se = math.sqrt(autocovariance(y, 0) / len(y))
    assert newey_west_mean_se(y, 30) > iid_se


def test_validation():
    x = [1.0, 2.0, 3.0, 4.0]
    with pytest.raises(ValueError):
        autocovariance(x, -1)
    with pytest.raises(ValueError):
        autocovariance(x, 4)                 # lag >= n
    with pytest.raises(ValueError):
        newey_west_variance(x, 4)            # lags >= n
    with pytest.raises(ValueError):
        autocorrelation([2.0, 2.0, 2.0], 1)  # zero variance
