"""Physical-measure (KMV) distance to default and default probability."""

import math
import random

import pytest

from quantforge import (physical_distance_to_default, physical_default_probability,
                        distance_to_default, risk_neutral_default_probability)


def test_mu_equals_r_recovers_risk_neutral():
    V, D, r, sigma, t = 120.0, 100.0, 0.03, 0.25, 1.0
    assert abs(physical_distance_to_default(V, D, r, sigma, t)
               - distance_to_default(V, D, r, sigma, t)) < 1e-12
    assert abs(physical_default_probability(V, D, r, sigma, t)
               - risk_neutral_default_probability(V, D, r, sigma, t)) < 1e-12


def test_physical_pd_below_risk_neutral_when_mu_above_r():
    V, D, r, sigma, t = 120.0, 100.0, 0.03, 0.25, 1.0
    mu = 0.10
    assert (physical_default_probability(V, D, mu, sigma, t)
            < risk_neutral_default_probability(V, D, r, sigma, t))


def test_matches_gbm_monte_carlo():
    V, D, mu, sigma, t = 120.0, 100.0, 0.10, 0.25, 1.0
    rng = random.Random(999)
    n = 500000
    d = 0
    for _ in range(n):
        z = rng.gauss(0, 1)
        vt = V * math.exp((mu - 0.5 * sigma * sigma) * t + sigma * math.sqrt(t) * z)
        if vt < D:
            d += 1
    assert abs(d / n - physical_default_probability(V, D, mu, sigma, t)) < 3e-3


def test_higher_drift_lowers_default_probability():
    V, D, sigma, t = 120.0, 100.0, 0.25, 1.0
    lo = physical_default_probability(V, D, 0.05, sigma, t)
    hi = physical_default_probability(V, D, 0.15, sigma, t)
    assert hi < lo


def test_pd_is_phi_of_negative_dd():
    from quantforge.mathfns import norm_cdf
    V, D, mu, sigma, t = 130.0, 90.0, 0.08, 0.3, 2.0
    dd = physical_distance_to_default(V, D, mu, sigma, t)
    assert abs(physical_default_probability(V, D, mu, sigma, t) - norm_cdf(-dd)) < 1e-14


def test_more_leverage_raises_default():
    mu, sigma, t = 0.08, 0.25, 1.0
    safe = physical_default_probability(150.0, 100.0, mu, sigma, t)
    risky = physical_default_probability(105.0, 100.0, mu, sigma, t)
    assert risky > safe


def test_validation():
    with pytest.raises(ValueError):
        physical_distance_to_default(0.0, 100.0, 0.08, 0.25, 1.0)
    with pytest.raises(ValueError):
        physical_default_probability(120.0, 100.0, 0.08, 0.0, 1.0)
    with pytest.raises(ValueError):
        physical_default_probability(120.0, 100.0, 0.08, 0.25, 0.0)
