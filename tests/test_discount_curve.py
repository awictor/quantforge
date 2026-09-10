"""Discount curve interpolation and par-swap bootstrap."""

import math

import pytest

from quantforge import DiscountCurve, bootstrap_from_swaps
from quantforge.g2pp import zero_bond as g2pp_zero_bond


def test_flat_curve_from_zero_rates():
    c = DiscountCurve.from_zero_rates([1, 2, 5, 10], [0.03] * 4)
    for T in (1, 3, 7, 10):
        assert c.df(T) == pytest.approx(math.exp(-0.03 * T), abs=1e-9)


def test_interpolation_exact_at_pillars():
    c = DiscountCurve([1, 3, 5], [0.97, 0.90, 0.83])
    assert c.df(1) == pytest.approx(0.97)
    assert c.df(3) == pytest.approx(0.90)
    assert c.df(5) == pytest.approx(0.83)


def test_df_monotone_decreasing_for_positive_rates():
    c = DiscountCurve.from_zero_rates([1, 5, 10], [0.03, 0.035, 0.04])
    prev = 1.0
    for T in [0.5, 1, 2, 5, 8, 10]:
        d = c.df(T)
        assert d < prev
        prev = d


def test_zero_and_forward_rate_consistency():
    c = DiscountCurve.from_zero_rates([1, 5, 10], [0.03, 0.035, 0.04])
    # Forward over [0, T] equals the zero rate.
    assert c.forward_rate(1e-9, 5) == pytest.approx(c.zero_rate(5), abs=1e-4)


def test_consecutive_tenor_bootstrap_is_exact():
    mats = [1, 2, 3, 4, 5, 6, 7]
    pars = [0.030, 0.032, 0.034, 0.035, 0.036, 0.0365, 0.037]
    curve = bootstrap_from_swaps(mats, pars, freq=1.0)
    for m, p in zip(mats, pars):
        pay = [k + 1 for k in range(m)]
        assert curve.par_swap_rate(pay) == pytest.approx(p, abs=1e-10)


def test_g2pp_reprices_curve_bonds_exactly():
    mats = [1, 2, 3, 4, 5]
    pars = [0.030, 0.032, 0.034, 0.035, 0.036]
    curve = bootstrap_from_swaps(mats, pars)
    a, b, sig, eta, rho = 0.1, 0.3, 0.01, 0.008, -0.3
    for T in (1, 3, 5):
        P = g2pp_zero_bond(curve.df(T), curve.df(0), 0.0, 0.0,
                           a, b, sig, eta, rho, 0.0, T)
        assert P == pytest.approx(curve.df(T), abs=1e-12)


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        DiscountCurve([1, 2], [0.97, -0.1])     # negative DF
    with pytest.raises(ValueError):
        DiscountCurve([-1, 2], [0.97, 0.9])     # negative maturity
    with pytest.raises(ValueError):
        bootstrap_from_swaps([1, 2], [0.03])    # mismatched lengths
