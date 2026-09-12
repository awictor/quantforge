"""Log contract and its variance-swap replication link."""

import math

import pytest

from quantforge import log_contract, log_contract_fair_variance


S, T, R, SIG = 100.0, 1.0, 0.05, 0.25


def test_value_negative():
    assert log_contract(S, T, R, SIG) < 0


def test_value_formula():
    assert log_contract(S, T, R, SIG) == pytest.approx(
        math.exp(-R * T) * (-0.5 * SIG ** 2 * T))


def test_fair_variance_recovers_sigma_squared():
    assert log_contract_fair_variance(S, T, R, SIG) == pytest.approx(SIG ** 2, abs=1e-12)


def test_more_negative_with_horizon_and_vol():
    assert log_contract(S, 2.0, R, SIG) < log_contract(S, 1.0, R, SIG)
    assert log_contract(S, T, R, 0.4) < log_contract(S, T, R, 0.25)


def test_zero_vol_zero_value():
    assert log_contract(S, T, R, 0.0) == pytest.approx(0.0, abs=1e-15)


def test_validation():
    with pytest.raises(ValueError):
        log_contract_fair_variance(S, 0, R, SIG)
