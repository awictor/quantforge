"""Inflation-linked bond and breakeven-inflation analytics."""

import pytest

from quantforge import (
    index_ratio, inflation_adjusted_principal, fisher_real_rate,
    fisher_nominal_rate, breakeven_inflation, real_from_breakeven,
    linker_price, linker_real_yield,
    bond_cashflows, bond_price_from_yield, yield_to_maturity,
)


CF = bond_cashflows(100, 0.03, 5, 2)  # real cashflows


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


def test_linker_price_ratio_one_matches_bond():
    # At ratio 1 the linker price is a plain real-yield bond price.
    assert linker_price(CF, 0.015, 100, 100) == pytest.approx(
        bond_price_from_yield(CF, 0.015), abs=1e-12)


def test_linker_price_homogeneous_in_ratio():
    base = bond_price_from_yield(CF, 0.015)
    assert linker_price(CF, 0.015, 120, 100) == pytest.approx(1.2 * base, abs=1e-10)


def test_linker_real_yield_inverts_price():
    p = linker_price(CF, 0.015, 120, 100)
    assert linker_real_yield(CF, p, 120, 100) == pytest.approx(0.015, abs=1e-8)


def test_linker_real_yield_equals_deflated_bond_ytm():
    p = linker_price(CF, 0.015, 120, 100)
    base = bond_price_from_yield(CF, 0.015)
    assert linker_real_yield(CF, p, 120, 100) == pytest.approx(
        yield_to_maturity(CF, base), abs=1e-8)


def test_linker_validation():
    with pytest.raises(ValueError):
        linker_price(CF, 0.015, 100, 0)
    with pytest.raises(ValueError):
        linker_real_yield(CF, -5, 100, 100)


def test_validation():
    with pytest.raises(ValueError):
        index_ratio(100, 0)
    with pytest.raises(ValueError):
        index_ratio(0, 100)
    with pytest.raises(ValueError):
        breakeven_inflation(0.05, -1.5)
    with pytest.raises(ValueError):
        real_from_breakeven(0.05, -1.5)
