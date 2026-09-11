"""Risk-neutral VaR and CVaR of the terminal return from a smile."""

import math
import statistics

import pytest

from quantforge import (
    risk_neutral_var_from_smile, risk_neutral_cvar_from_smile,
)
from quantforge.mathfns import norm_cdf


S0, T, R, SIG = 100.0, 1.0, 0.03, 0.2
_N = statistics.NormalDist()


def _flat(K):
    return SIG


def _lognormal_var_cvar(alpha):
    m = math.log(S0) + (R - 0.5 * SIG * SIG) * T
    s = SIG * math.sqrt(T)
    tail = 1.0 - alpha
    K = math.exp(m + s * _N.inv_cdf(tail))
    var = 1.0 - K / S0
    es = math.exp(m + 0.5 * s * s) * norm_cdf((math.log(K) - m - s * s) / s)
    cvar = 1.0 - es / (S0 * tail)
    return var, cvar


@pytest.mark.parametrize("alpha", [0.90, 0.95, 0.99])
def test_var_matches_lognormal(alpha):
    var_cf, _ = _lognormal_var_cvar(alpha)
    assert risk_neutral_var_from_smile(S0, T, R, _flat, alpha) == pytest.approx(
        var_cf, abs=1e-3)


@pytest.mark.parametrize("alpha", [0.90, 0.95, 0.99])
def test_cvar_matches_lognormal(alpha):
    _, cvar_cf = _lognormal_var_cvar(alpha)
    assert risk_neutral_cvar_from_smile(S0, T, R, _flat, alpha) == pytest.approx(
        cvar_cf, abs=1e-3)


def test_cvar_exceeds_var():
    for alpha in (0.90, 0.95, 0.99):
        var = risk_neutral_var_from_smile(S0, T, R, _flat, alpha)
        cvar = risk_neutral_cvar_from_smile(S0, T, R, _flat, alpha)
        assert cvar >= var


def test_var_increases_with_confidence():
    v90 = risk_neutral_var_from_smile(S0, T, R, _flat, 0.90)
    v99 = risk_neutral_var_from_smile(S0, T, R, _flat, 0.99)
    assert v99 > v90


def test_skew_raises_tail_risk():
    # Equity skew (higher vol at low strikes) fattens the loss tail, so both VaR
    # and CVaR rise versus a flat smile at the same ATM vol.
    def skew(K):
        return max(0.05, SIG + 0.2 * math.log(S0 / K))

    for alpha in (0.95, 0.99):
        assert (risk_neutral_var_from_smile(S0, T, R, skew, alpha)
                > risk_neutral_var_from_smile(S0, T, R, _flat, alpha))
        assert (risk_neutral_cvar_from_smile(S0, T, R, skew, alpha)
                > risk_neutral_cvar_from_smile(S0, T, R, _flat, alpha))


def test_alpha_out_of_range_raises():
    with pytest.raises(ValueError):
        risk_neutral_var_from_smile(S0, T, R, _flat, 1.0)
    with pytest.raises(ValueError):
        risk_neutral_cvar_from_smile(S0, T, R, _flat, 0.0)
