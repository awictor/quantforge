"""Tests for CEV pricing and the noncentral chi-square CDF."""

import math

import pytest

from quantforge import cev_price, noncentral_chisq_cdf, call_price, put_price, OptionType
from quantforge.cev import _gammainc_lower_reg


# --- noncentral chi-square CDF ---
def test_gammainc_known_values():
    # P(1, x) = 1 - e^{-x}.
    assert _gammainc_lower_reg(1.0, 2.0) == pytest.approx(1 - math.exp(-2), abs=1e-12)


def test_noncentral_reduces_to_central_at_zero_lambda():
    # lam = 0 -> central chi-square CDF = P(k/2, x/2).
    assert noncentral_chisq_cdf(3.0, 4.0, 0.0) == pytest.approx(
        _gammainc_lower_reg(2.0, 1.5), abs=1e-12)


def test_noncentral_stable_at_large_lambda():
    # At the distribution mean (= k + lam) the CDF is近 0.5, even for huge lam
    # where a naive j=0 Poisson start would underflow.
    assert noncentral_chisq_cdf(5010, 10, 5000) == pytest.approx(0.5, abs=0.02)


def test_noncentral_monotone_in_x():
    prev = -1.0
    for x in (1, 5, 10, 20, 40):
        v = noncentral_chisq_cdf(x, 4, 3)
        assert v >= prev
        prev = v


# --- CEV pricing ---
@pytest.mark.parametrize("beta", [0.5, 0.7, 0.9, 0.99])
def test_cev_matches_bsm_when_scaled(beta):
    # delta is calibrated so the ATM level matches sigma, so CEV stays close to
    # Black-Scholes across beta for an ATM option.
    c = cev_price(100, 100, 1.0, 0.05, 0.2, beta, OptionType.CALL)
    assert c == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2), abs=0.05)


def test_put_call_parity():
    c = cev_price(100, 95, 1.0, 0.05, 0.25, 0.6, OptionType.CALL)
    p = cev_price(100, 95, 1.0, 0.05, 0.25, 0.6, OptionType.PUT)
    assert c - p == pytest.approx(100 - 95 * math.exp(-0.05), abs=1e-4)


def test_lower_beta_richer_otm_put():
    # Lower beta => stronger leverage/skew => a downside put is worth more.
    low = cev_price(100, 85, 1.0, 0.05, 0.25, 0.3, OptionType.PUT)
    high = cev_price(100, 85, 1.0, 0.05, 0.25, 0.9, OptionType.PUT)
    assert low > high


def test_price_positive_and_bounded():
    c = cev_price(100, 100, 0.5, 0.03, 0.3, 0.5, OptionType.CALL)
    assert 0 < c < 100


def test_zero_time_is_intrinsic():
    assert cev_price(110, 100, 0.0, 0.05, 0.2, 0.5, OptionType.CALL) == pytest.approx(10.0)


def test_rejects_beta_out_of_range():
    with pytest.raises(ValueError):
        cev_price(100, 100, 1.0, 0.05, 0.2, 1.0)   # beta must be < 1
    with pytest.raises(ValueError):
        cev_price(100, 100, 1.0, 0.05, 0.2, -0.1)
