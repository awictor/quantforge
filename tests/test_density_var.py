"""Option-position VaR/ES under the smile-implied risk-neutral density."""

import math

import pytest

from quantforge import density_var_es, call_price
from quantforge.mathfns import norm_ppf


S0, T, R = 100.0, 1.0, 0.05
SIG = 0.2
F = S0 * math.exp(R * T)


def _flat(K):
    return SIG


def test_long_forward_var_matches_lognormal_quantile():
    z = norm_ppf(0.01)
    sT = F * math.exp(-0.5 * SIG * SIG * T + SIG * math.sqrt(T) * z)
    var_analytic = F - sT
    var, _es = density_var_es(S0, T, R, _flat, lambda ST: ST - F,
                              confidence=0.99)
    assert var == pytest.approx(var_analytic, abs=1.0)


def test_es_at_least_var():
    var, es = density_var_es(S0, T, R, _flat, lambda ST: ST - F,
                             confidence=0.99)
    assert es >= var


def test_long_call_var_capped_at_premium():
    prem = call_price(S0, 100, T, R, SIG)
    var, es = density_var_es(S0, T, R, _flat, lambda ST: max(ST - 100, 0) - prem,
                             confidence=0.99)
    # Long option's max loss is the premium.
    assert var == pytest.approx(prem, abs=1e-1)
    assert es <= prem + 1e-6


def test_short_call_var_positive_and_es_larger():
    prem = call_price(S0, 100, T, R, SIG)
    var, es = density_var_es(S0, T, R, _flat, lambda ST: prem - max(ST - 100, 0),
                             confidence=0.95)
    assert var > 0.0
    assert es >= var


def test_higher_confidence_larger_var():
    lo, _ = density_var_es(S0, T, R, _flat, lambda ST: ST - F, confidence=0.95)
    hi, _ = density_var_es(S0, T, R, _flat, lambda ST: ST - F, confidence=0.99)
    assert hi > lo
