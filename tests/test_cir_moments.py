"""CIR short-rate moments: expected rate, variance, stationary Gamma law."""

import math

import pytest

from quantforge import (
    cir_expected_rate, cir_rate_variance, cir_stationary_distribution,
)


R0, KAPPA, THETA, SIG = 0.03, 0.5, 0.05, 0.05


def test_expected_rate_formula():
    t = 2.0
    exp = THETA + (R0 - THETA) * math.exp(-KAPPA * t)
    assert cir_expected_rate(R0, t, KAPPA, THETA) == pytest.approx(exp, abs=1e-12)


@pytest.mark.slow
def test_moments_match_monte_carlo():
    import random
    rng = random.Random(3)
    t, N, nsteps = 2.0, 60000, 400
    dt = t / nsteps
    vals = []
    for _ in range(N):
        r = R0
        for _ in range(nsteps):
            r += KAPPA * (THETA - r) * dt + SIG * math.sqrt(max(r, 0.0) * dt) * rng.gauss(0, 1)
            r = max(r, 0.0)
        vals.append(r)
    mean = sum(vals) / N
    var = sum((v - mean) ** 2 for v in vals) / N
    assert cir_expected_rate(R0, t, KAPPA, THETA) == pytest.approx(mean, abs=1e-3)
    assert cir_rate_variance(R0, t, KAPPA, THETA, SIG) == pytest.approx(var, rel=0.05)


def test_variance_approaches_stationary():
    _, _, _, svar = cir_stationary_distribution(KAPPA, THETA, SIG)
    assert cir_rate_variance(R0, 0.0, KAPPA, THETA, SIG) == pytest.approx(0.0, abs=1e-15)
    assert cir_rate_variance(R0, 500.0, KAPPA, THETA, SIG) == pytest.approx(svar, abs=1e-9)


def test_stationary_gamma_consistency():
    shape, scale, mean, var = cir_stationary_distribution(KAPPA, THETA, SIG)
    assert shape * scale == pytest.approx(mean, abs=1e-12)
    assert shape * scale * scale == pytest.approx(var, abs=1e-12)
    assert mean == pytest.approx(THETA, abs=1e-12)


def test_variance_depends_on_r0():
    # Unlike Vasicek, the CIR variance rises with the starting rate.
    lo = cir_rate_variance(0.01, 1.0, KAPPA, THETA, SIG)
    hi = cir_rate_variance(0.10, 1.0, KAPPA, THETA, SIG)
    assert hi > lo


def test_validation():
    with pytest.raises(ValueError):
        cir_expected_rate(R0, 1.0, 0.0, THETA)
    with pytest.raises(ValueError):
        cir_stationary_distribution(KAPPA, THETA, 0.0)
