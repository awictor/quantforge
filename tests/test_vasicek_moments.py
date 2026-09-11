"""Vasicek short-rate moments: expected rate, variance, stationary law."""

import math

import pytest

from quantforge import (
    vasicek_expected_rate, vasicek_rate_variance,
    vasicek_stationary_distribution,
)


R0, KAPPA, THETA, SIG = 0.03, 0.5, 0.05, 0.01


def test_expected_rate_formula():
    t = 2.0
    exp = THETA + (R0 - THETA) * math.exp(-KAPPA * t)
    assert vasicek_expected_rate(R0, t, KAPPA, THETA) == pytest.approx(exp, abs=1e-12)


@pytest.mark.slow
def test_expected_rate_matches_monte_carlo():
    import random
    rng = random.Random(3)
    t, N, nsteps = 2.0, 120000, 300
    dt = t / nsteps
    total = 0.0
    for _ in range(N):
        r = R0
        for _ in range(nsteps):
            r += KAPPA * (THETA - r) * dt + SIG * math.sqrt(dt) * rng.gauss(0, 1)
        total += r
    mc = total / N
    assert vasicek_expected_rate(R0, t, KAPPA, THETA) == pytest.approx(mc, abs=1e-3)


def test_variance_grows_to_stationary():
    _, svar = vasicek_stationary_distribution(KAPPA, THETA, SIG)
    assert vasicek_rate_variance(0.0, KAPPA, SIG) == pytest.approx(0.0, abs=1e-15)
    assert vasicek_rate_variance(1000.0, KAPPA, SIG) == pytest.approx(svar, abs=1e-9)
    assert svar == pytest.approx(SIG * SIG / (2 * KAPPA), abs=1e-12)


def test_zero_kappa_limits():
    t = 2.0
    assert vasicek_expected_rate(R0, t, 0.0, THETA) == pytest.approx(R0, abs=1e-12)
    assert vasicek_rate_variance(t, 0.0, SIG) == pytest.approx(SIG * SIG * t, abs=1e-12)


def test_stationary_mean_is_theta():
    mean, _ = vasicek_stationary_distribution(KAPPA, THETA, SIG)
    assert mean == pytest.approx(THETA, abs=1e-12)


def test_stationary_requires_positive_kappa():
    with pytest.raises(ValueError):
        vasicek_stationary_distribution(0.0, THETA, SIG)
