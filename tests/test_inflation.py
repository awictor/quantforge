"""Inflation-linked bond and breakeven-inflation analytics."""

import pytest

from quantforge import (
    index_ratio, inflation_adjusted_principal, fisher_real_rate,
    fisher_nominal_rate, breakeven_inflation, real_from_breakeven,
)


def test_index_ratio():
    assert index_ratio(110, 100) == pytest.approx(1.1)
    assert index_ratio(100, 100) == pytest.approx(1.0)


def test_inflation_adjusted_principal():
    assert inflation_adjusted_principal(1000, 110, 100) == pytest.approx(1100.0)


def test_fisher_round_trip():
    # real -> nominal -> real is identity for fixed inflation.
    real = fisher_real_rate(0.05, 0.02)
    assert fisher_nominal_rate(real, 0.02) == pytest.approx(0.05, abs=1e-12)


def test_fisher_real_below_nominal_minus_inflation():
    # Exact real rate is slightly below the linear approx nominal - inflation.
    real = fisher_real_rate(0.05, 0.02)
    assert real < 0.03
    assert real == pytest.approx(0.03 / 1.02, abs=1e-12)


def test_breakeven_recovers_inflation():
    # nominal from Fisher(real=2%, infl=3%) -> breakeven vs 2% real == 3%.
    nominal = fisher_nominal_rate(0.02, 0.03)
    assert breakeven_inflation(nominal, 0.02) == pytest.approx(0.03, abs=1e-12)


def test_breakeven_real_inverse():
    nominal = 0.05
    bkv = breakeven_inflation(nominal, 0.02)
    assert real_from_breakeven(nominal, bkv) == pytest.approx(0.02, abs=1e-12)


def test_breakeven_approx_difference_small_rates():
    # For small rates breakeven ~ nominal - real.
    bkv = breakeven_inflation(0.041, 0.02)
    assert bkv == pytest.approx(0.021, abs=5e-4)


def test_validation():
    with pytest.raises(ValueError):
        index_ratio(100, 0)
    with pytest.raises(ValueError):
        index_ratio(0, 100)
    with pytest.raises(ValueError):
        breakeven_inflation(0.05, -1.5)
    with pytest.raises(ValueError):
        real_from_breakeven(0.05, -1.5)
