"""Empirical historical VaR / CVaR from a return series (perfmetrics module)."""

import random

import pytest

from quantforge import historical_var_series, historical_cvar


def _series(n=5000, seed=3):
    rng = random.Random(seed)
    return [rng.gauss(0.0005, 0.01) for _ in range(n)]


def test_positive_loss():
    assert historical_var_series(_series()) > 0.0


def test_cvar_at_least_var():
    r = _series()
    assert historical_cvar(r) >= historical_var_series(r)


def test_higher_confidence_larger():
    r = _series()
    assert historical_var_series(r, 0.99) > historical_var_series(r, 0.95)
    assert historical_cvar(r, 0.99) > historical_cvar(r, 0.95)


def test_tail_fraction_matches_confidence():
    r = _series()
    v = historical_var_series(r, 0.95)
    frac = sum(1 for x in r if x <= -v) / len(r)
    assert frac == pytest.approx(0.05, abs=0.01)


def test_matches_known_percentile():
    # Uniform-ish small sample: 95% VaR is the negated 5th percentile.
    r = [-0.10, -0.05, -0.02, 0.0, 0.01, 0.03, 0.05, 0.08, 0.10, 0.12]
    # historical_var_series uses linear-interpolated 5th percentile.
    assert historical_var_series(r, 0.95) > 0.0


def test_validation():
    with pytest.raises(ValueError):
        historical_var_series([0.01])
    with pytest.raises(ValueError):
        historical_cvar([0.01])
