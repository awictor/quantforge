"""Sample risk measures: VaR, ES, spectral, entropic."""

import pytest

from quantforge import (
    value_at_risk, sample_expected_shortfall, spectral_risk_exponential,
    entropic_risk,
)


PNL = [0.01, -0.02, 0.03, -0.05, 0.02, -0.10, 0.015, -0.03, 0.025, -0.08,
       0.01, -0.04, 0.02, -0.06, 0.005, -0.12, 0.03, -0.02, 0.01, -0.07]


def test_es_at_least_var():
    assert sample_expected_shortfall(PNL, 0.95) >= value_at_risk(PNL, 0.95)


def test_var_increasing_in_confidence():
    assert value_at_risk(PNL, 0.99) >= value_at_risk(PNL, 0.90)


def test_spectral_increasing_in_risk_aversion():
    assert spectral_risk_exponential(PNL, 10.0) > spectral_risk_exponential(PNL, 2.0)


def test_entropic_approaches_mean_loss():
    mean_loss = -sum(PNL) / len(PNL)
    assert entropic_risk(PNL, 0.001) == pytest.approx(mean_loss, abs=1e-3)


def test_entropic_increasing_in_risk_aversion():
    assert entropic_risk(PNL, 5.0) > entropic_risk(PNL, 1.0)


def test_validation():
    with pytest.raises(ValueError):
        value_at_risk(PNL, 1.5)
    with pytest.raises(ValueError):
        sample_expected_shortfall([], 0.95)
    with pytest.raises(ValueError):
        entropic_risk(PNL, 0)
