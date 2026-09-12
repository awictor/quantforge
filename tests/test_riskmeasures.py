"""Sample risk measures: VaR, ES, spectral, entropic."""

import pytest

import random

from quantforge import (
    value_at_risk, sample_expected_shortfall, spectral_risk_exponential,
    entropic_risk, is_subadditive, component_expected_shortfall,
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


def test_expected_shortfall_subadditive():
    random.seed(2)
    a = [random.gauss(0, 0.02) for _ in range(200)]
    b = [random.gauss(0, 0.03) for _ in range(200)]
    assert is_subadditive(a, b, 0.95)
    assert is_subadditive(a, a, 0.95)   # perfectly correlated: still <=


def test_component_es_sums_to_total():
    random.seed(2)
    a = [random.gauss(0, 0.02) for _ in range(200)]
    b = [random.gauss(0, 0.03) for _ in range(200)]
    total = [a[i] + b[i] for i in range(200)]
    contribs = component_expected_shortfall([a, b], 0.95)
    assert sum(contribs) == pytest.approx(sample_expected_shortfall(total, 0.95),
                                          abs=1e-9)


def test_component_es_three_components():
    random.seed(3)
    a = [random.gauss(0, 0.02) for _ in range(200)]
    b = [random.gauss(0, 0.03) for _ in range(200)]
    c = [random.gauss(0, 0.01) for _ in range(200)]
    total = [a[i] + b[i] + c[i] for i in range(200)]
    assert sum(component_expected_shortfall([a, b, c], 0.95)) == pytest.approx(
        sample_expected_shortfall(total, 0.95), abs=1e-9)


def test_coherence_validation():
    with pytest.raises(ValueError):
        is_subadditive([0.1, 0.2], [0.1])
    with pytest.raises(ValueError):
        component_expected_shortfall([])


def test_validation():
    with pytest.raises(ValueError):
        value_at_risk(PNL, 1.5)
    with pytest.raises(ValueError):
        sample_expected_shortfall([], 0.95)
    with pytest.raises(ValueError):
        entropic_risk(PNL, 0)
