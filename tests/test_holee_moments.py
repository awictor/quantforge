"""Ho-Lee short-rate moments: expected rate and variance (drifted Brownian)."""

import math

import pytest

from quantforge import holee_expected_rate, holee_rate_variance


R0, THETA, SIG = 0.03, 0.01, 0.02


def test_expected_rate_is_linear_drift():
    for t in (0.5, 1.0, 3.0):
        assert holee_expected_rate(R0, t, THETA) == pytest.approx(R0 + THETA * t, abs=1e-12)


def test_variance_is_linear_in_time():
    for t in (0.5, 1.0, 3.0):
        assert holee_rate_variance(t, SIG) == pytest.approx(SIG * SIG * t, abs=1e-12)


@pytest.mark.slow
def test_moments_match_monte_carlo():
    import random
    rng = random.Random(3)
    t, N, nsteps = 3.0, 60000, 300
    dt = t / nsteps
    vals = []
    for _ in range(N):
        r = R0
        for _ in range(nsteps):
            r += THETA * dt + SIG * math.sqrt(dt) * rng.gauss(0, 1)
        vals.append(r)
    mean = sum(vals) / N
    var = sum((v - mean) ** 2 for v in vals) / N
    assert holee_expected_rate(R0, t, THETA) == pytest.approx(mean, abs=1e-3)
    assert holee_rate_variance(t, SIG) == pytest.approx(var, rel=0.03)


def test_variance_unbounded():
    # No mean reversion: variance keeps growing (no stationary cap).
    assert holee_rate_variance(100.0, SIG) > holee_rate_variance(10.0, SIG)


def test_validation():
    with pytest.raises(ValueError):
        holee_expected_rate(R0, -1.0, THETA)
    with pytest.raises(ValueError):
        holee_rate_variance(-1.0, SIG)
