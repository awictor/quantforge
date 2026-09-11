"""Inflation-linked bond and breakeven-inflation analytics."""

import pytest

from quantforge import (
    index_ratio, inflation_adjusted_principal, fisher_real_rate,
    fisher_nominal_rate, breakeven_inflation, real_from_breakeven,
    linker_price, linker_real_yield,
    linker_real_duration, linker_real_convexity, linker_real_dv01,
    deflation_floored_redemption, deflation_floor_value, yoy_inflation_rate,
    zc_inflation_swap_rate, zc_inflation_swap_value,
    inflation_curve_from_zc_swaps, forward_inflation_rate, yoy_swap_value,
    reference_cpi, index_ratio_interpolated,
    normalize_seasonal_factors, apply_seasonality, deseasonalize,
    bond_cashflows, bond_price_from_yield, yield_to_maturity,
)


RAW_SEASONAL = [1.02, 0.99, 1.01, 1.00, 0.98, 1.03,
                1.01, 0.97, 1.02, 1.00, 0.99, 1.01]
from quantforge import modified_duration, convexity, bond_dv01


TENORS = [1, 2, 3, 5]
ZC = [0.025, 0.028, 0.03, 0.032]


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


def test_floor_ootm_in_inflation():
    # Ratio > 1: redemption equals adjusted principal, floor worthless.
    assert deflation_floored_redemption(1000, 120, 100) == pytest.approx(
        inflation_adjusted_principal(1000, 120, 100))
    assert deflation_floor_value(1000, 120, 100) == pytest.approx(0.0, abs=1e-12)


def test_floor_binds_in_deflation():
    assert deflation_floored_redemption(1000, 90, 100) == pytest.approx(1000.0)
    assert deflation_floor_value(1000, 90, 100) == pytest.approx(100.0)


def test_floor_never_below_par():
    for idx in (60, 90, 100, 150):
        assert deflation_floored_redemption(1000, idx, 100) >= 1000 - 1e-9


def test_yoy():
    assert yoy_inflation_rate(100, 103) == pytest.approx(0.03)


def test_zc_swap_rate_compounding_identity():
    k = zc_inflation_swap_rate(100, 133.1, 3)
    assert (1 + k) ** 3 * 100 == pytest.approx(133.1, abs=1e-9)


def test_zc_swap_value_zero_at_par():
    k = zc_inflation_swap_rate(100, 133.1, 3)
    assert zc_inflation_swap_value(1e6, k, 100, 133.1, 3) == pytest.approx(0.0, abs=1e-3)


def test_zc_swap_receiver_gains_above_par():
    assert zc_inflation_swap_value(1e6, 0.05, 100, 133.1, 3) > 0


def test_inflation_swap_validation():
    with pytest.raises(ValueError):
        zc_inflation_swap_rate(100, 133, 0)
    with pytest.raises(ValueError):
        yoy_inflation_rate(0, 100)


def test_curve_reprices_input_swaps():
    lv = inflation_curve_from_zc_swaps(100.0, TENORS, ZC)
    for i, T in enumerate(TENORS):
        assert zc_inflation_swap_rate(100.0, lv[i], T) == pytest.approx(ZC[i], abs=1e-12)


def test_forward_chains_to_spot():
    lv = inflation_curve_from_zc_swaps(100.0, TENORS, ZC)
    f12 = forward_inflation_rate(lv[0], lv[1], 1, 2)
    assert (1 + ZC[0]) ** 1 * (1 + f12) ** 1 == pytest.approx((1 + ZC[1]) ** 2, abs=1e-9)


def test_forward_mid_curve_chain():
    lv = inflation_curve_from_zc_swaps(100.0, TENORS, ZC)
    f23 = forward_inflation_rate(lv[1], lv[2], 2, 3)
    assert (1 + ZC[1]) ** 2 * (1 + f23) == pytest.approx((1 + ZC[2]) ** 3, abs=1e-9)


def test_yoy_swap_zero_when_fixed_equals_flat_inflation():
    lv = [100 * 1.03 ** t for t in (1, 2, 3)]
    assert yoy_swap_value(1e6, 0.03, lv, [1, 1, 1], 100) == pytest.approx(0.0, abs=1e-6)


def test_yoy_swap_receiver_gains_below_fixed():
    lv = [100 * 1.03 ** t for t in (1, 2, 3)]
    assert yoy_swap_value(1e6, 0.02, lv, [1, 1, 1], 100) > 0


def test_curve_forward_validation():
    with pytest.raises(ValueError):
        forward_inflation_rate(100, 110, 2, 2)
    with pytest.raises(ValueError):
        inflation_curve_from_zc_swaps(100, [1, 2], [0.02])
    with pytest.raises(ValueError):
        yoy_swap_value(1e6, 0.02, [110], [1, 1], 100)


def test_reference_cpi_hits_month_start():
    assert reference_cpi(200, 203, 1, 30) == pytest.approx(200.0)


def test_reference_cpi_linear_midpoint():
    # Day 16 of 30 -> frac 15/30 = 0.5 -> mean of anchors.
    assert reference_cpi(200, 203, 16, 30) == pytest.approx(201.5)


def test_reference_cpi_monotone():
    vals = [reference_cpi(200, 203, d, 30) for d in range(1, 31)]
    assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))


def test_interpolated_ratio_matches_ref_over_base():
    assert index_ratio_interpolated(200, 203, 16, 30, 190) == pytest.approx(201.5 / 190)


def test_interpolated_ratio_day1_equals_plain():
    assert index_ratio_interpolated(200, 203, 1, 30, 190) == pytest.approx(
        index_ratio(200, 190))


def test_reference_cpi_validation():
    with pytest.raises(ValueError):
        reference_cpi(200, 203, 31, 30)
    with pytest.raises(ValueError):
        reference_cpi(200, 203, 0, 30)
    with pytest.raises(ValueError):
        reference_cpi(-1, 203, 1, 30)


def test_real_duration_matches_bondmath():
    assert linker_real_duration(CF, 0.015) == pytest.approx(
        modified_duration(CF, 0.015), abs=1e-12)


def test_real_convexity_matches_bondmath():
    assert linker_real_convexity(CF, 0.015) == pytest.approx(
        convexity(CF, 0.015), abs=1e-12)


def test_real_duration_index_independent():
    # Fractional duration cancels the index ratio.
    d1 = linker_real_duration(CF, 0.015)
    assert d1 == modified_duration(CF, 0.015)


def test_real_duration_finite_difference():
    r, h = 0.015, 1e-6
    P = lambda y: linker_price(CF, y, 120, 100)
    fd = -(P(r + h) - P(r - h)) / (2 * h) / P(r)
    assert linker_real_duration(CF, r) == pytest.approx(fd, abs=1e-5)


def test_real_convexity_finite_difference():
    r, h = 0.015, 1e-6
    P = lambda y: linker_price(CF, y, 120, 100)
    fd = (P(r + h) - 2 * P(r) + P(r - h)) / h ** 2 / P(r)
    assert linker_real_convexity(CF, r) == pytest.approx(fd, abs=1e-3)


def test_real_dv01_scales_with_ratio():
    assert linker_real_dv01(CF, 0.015, 120, 100) == pytest.approx(
        1.2 * bond_dv01(CF, 0.015), abs=1e-9)
    assert linker_real_dv01(CF, 0.015, 120, 100) < 0


def test_real_risk_validation():
    with pytest.raises(ValueError):
        linker_real_duration([(1, -100)], 5.0)


def test_seasonal_factors_product_one():
    f = normalize_seasonal_factors(RAW_SEASONAL)
    prod = 1.0
    for x in f:
        prod *= x
    assert prod == pytest.approx(1.0, abs=1e-12)


def test_seasonal_normalization_idempotent():
    f = normalize_seasonal_factors(RAW_SEASONAL)
    f2 = normalize_seasonal_factors(f)
    assert all(f[i] == pytest.approx(f2[i], abs=1e-12) for i in range(12))


def test_seasonal_preserves_relative_shape():
    f = normalize_seasonal_factors(RAW_SEASONAL)
    assert f[0] / f[1] == pytest.approx(RAW_SEASONAL[0] / RAW_SEASONAL[1])


def test_apply_deseasonalize_roundtrip():
    f = normalize_seasonal_factors(RAW_SEASONAL)
    assert deseasonalize(apply_seasonality(250.0, f[3]), f[3]) == pytest.approx(250.0)


def test_seasonal_validation():
    with pytest.raises(ValueError):
        normalize_seasonal_factors([1, 2, 3])
    with pytest.raises(ValueError):
        normalize_seasonal_factors([1] * 11 + [-1])
    with pytest.raises(ValueError):
        apply_seasonality(-1, 1.0)
    with pytest.raises(ValueError):
        deseasonalize(100, 0)


def test_validation():
    with pytest.raises(ValueError):
        index_ratio(100, 0)
    with pytest.raises(ValueError):
        index_ratio(0, 100)
    with pytest.raises(ValueError):
        breakeven_inflation(0.05, -1.5)
    with pytest.raises(ValueError):
        real_from_breakeven(0.05, -1.5)
