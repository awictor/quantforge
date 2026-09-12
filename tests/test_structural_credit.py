"""Merton structural credit model."""

import pytest

from quantforge import (
    equity_value, risk_neutral_default_probability, distance_to_default,
    risky_debt_value, credit_spread,
)
from quantforge.bsm import call_price
from quantforge.mathfns import norm_cdf


V, D, R, SIG, T = 120.0, 100.0, 0.05, 0.25, 1.0


def test_equity_is_call_on_assets():
    assert equity_value(V, D, R, SIG, T) == pytest.approx(
        call_price(V, D, T, R, SIG, b=R), abs=1e-12)


def test_pd_is_phi_minus_distance():
    assert risk_neutral_default_probability(V, D, R, SIG, T) == pytest.approx(
        norm_cdf(-distance_to_default(V, D, R, SIG, T)), abs=1e-12)


def test_value_identity():
    assert V == pytest.approx(
        equity_value(V, D, R, SIG, T) + risky_debt_value(V, D, R, SIG, T), abs=1e-9)


def test_spread_positive():
    assert credit_spread(V, D, R, SIG, T) > 0


def test_leverage_raises_pd_and_spread():
    assert risk_neutral_default_probability(105, D, R, SIG, T) > \
        risk_neutral_default_probability(V, D, R, SIG, T)
    assert credit_spread(105, D, R, SIG, T) > credit_spread(V, D, R, SIG, T)


def test_volatility_raises_pd_and_spread():
    assert risk_neutral_default_probability(V, D, R, 0.4, T) > \
        risk_neutral_default_probability(V, D, R, SIG, T)
    assert credit_spread(V, D, R, 0.4, T) > credit_spread(V, D, R, SIG, T)


def test_safe_firm_tiny_risk():
    assert risk_neutral_default_probability(300, D, R, 0.1, T) < 1e-6
    assert credit_spread(300, D, R, 0.1, T) < 1e-4


def test_validation():
    with pytest.raises(ValueError):
        equity_value(-1, D, R, SIG, T)
    with pytest.raises(ValueError):
        distance_to_default(V, D, R, SIG, 0)
