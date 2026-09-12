"""Warrants and employee stock options."""

import math

import pytest

from quantforge import (
    dilution_factor, warrant_price, eso_expected_life, eso_value,
)
from quantforge.bsm import call_price


S, K, T, R, SIG = 50.0, 50.0, 5.0, 0.05, 0.3


def test_dilution_factor():
    assert dilution_factor(1_000_000, 100_000) == pytest.approx(1e6 / 1.1e6)
    assert dilution_factor(1e6, 0) == 1.0


def test_warrant_below_vanilla_call():
    van = call_price(S, K, T, R, SIG)
    w = warrant_price(S, K, T, R, SIG, 1_000_000, 100_000)
    assert w < van
    assert w / van == pytest.approx(dilution_factor(1_000_000, 100_000))


def test_expected_life_between_vesting_and_term():
    life = eso_expected_life(2.0, 10.0, 0.15)
    assert 2.0 < life < 10.0


def test_expected_life_no_exit_is_full_term():
    assert eso_expected_life(2.0, 10.0, 0.0) == 10.0


def test_eso_below_vanilla_on_full_term():
    eso = eso_value(S, K, 10.0, R, SIG, 2.0, 0.15, 0.03)
    assert eso < call_price(S, K, 10.0, R, SIG)


def test_eso_reduces_to_vanilla():
    assert eso_value(S, K, 10.0, R, SIG, 2.0, 0.0, 0.0) == pytest.approx(
        call_price(S, K, 10.0, R, SIG), abs=1e-9)


def test_higher_exit_rate_lowers_value():
    assert eso_value(S, K, 10.0, R, SIG, 2.0, 0.3, 0.0) < \
        eso_value(S, K, 10.0, R, SIG, 2.0, 0.1, 0.0)


def test_forfeiture_lowers_value():
    assert eso_value(S, K, 10.0, R, SIG, 2.0, 0.15, 0.1) < \
        eso_value(S, K, 10.0, R, SIG, 2.0, 0.15, 0.0)


def test_validation():
    with pytest.raises(ValueError):
        dilution_factor(-1, 100)
    with pytest.raises(ValueError):
        eso_expected_life(11, 10, 0.1)
    with pytest.raises(ValueError):
        eso_value(S, K, 10.0, R, SIG, 2.0, 0.15, -1)
