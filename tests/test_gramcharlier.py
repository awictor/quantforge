"""Tests for Corrado-Su skew/kurtosis pricing and realized moments."""

import math
import random

import pytest

from quantforge import (
    corrado_su_call, corrado_su_price, realized_skewness, realized_excess_kurtosis,
    call_price, put_price, OptionType,
)


def test_reduces_to_bsm_at_zero_moments():
    for K in (80, 100, 120):
        cs = corrado_su_call(100, K, 1.0, 0.05, 0.2, skew=0.0, excess_kurt=0.0)
        assert cs == pytest.approx(call_price(100, K, 1.0, 0.05, 0.2), abs=1e-9)


def test_put_call_parity():
    kw = dict(S=100, K=95, t=1.0, r=0.05, sigma=0.25)
    c = corrado_su_price(**kw, skew=-0.5, excess_kurt=2.0, option_type=OptionType.CALL)
    p = corrado_su_price(**kw, skew=-0.5, excess_kurt=2.0, option_type=OptionType.PUT)
    assert c - p == pytest.approx(100 - 95 * math.exp(-0.05), abs=1e-9)


def test_kurtosis_raises_far_wing_lowers_near_money():
    # Gram-Charlier kurtosis fattens the tails: it richens deep-wing options and
    # cheapens near-the-money ones. Use a short tenor so strikes sit far in
    # standardized-moneyness (d) space.
    t, r, sigma = 0.25, 0.0, 0.15
    far = 130  # deep OTM: d well below zero
    base_far = corrado_su_call(100, far, t, r, sigma, 0.0, 0.0)
    fat_far = corrado_su_call(100, far, t, r, sigma, 0.0, excess_kurt=3.0)
    assert fat_far > base_far

    base_atm = corrado_su_call(100, 100, t, r, sigma, 0.0, 0.0)
    fat_atm = corrado_su_call(100, 100, t, r, sigma, 0.0, excess_kurt=3.0)
    assert fat_atm < base_atm


def test_skew_shifts_call_value():
    # Nonzero skew moves the price away from Black-Scholes.
    base = corrado_su_call(100, 110, 1.0, 0.05, 0.2, 0.0, 0.0)
    skewed = corrado_su_call(100, 110, 1.0, 0.05, 0.2, skew=0.8, excess_kurt=0.0)
    assert skewed != pytest.approx(base)


def test_realized_skewness_of_normal_is_zero():
    rng = random.Random(1)
    xs = [rng.gauss(0, 1) for _ in range(100_000)]
    assert realized_skewness(xs) == pytest.approx(0.0, abs=0.05)


def test_realized_excess_kurtosis_of_normal_is_zero():
    rng = random.Random(2)
    xs = [rng.gauss(0, 1) for _ in range(100_000)]
    assert realized_excess_kurtosis(xs) == pytest.approx(0.0, abs=0.1)


def test_realized_excess_kurtosis_of_uniform():
    rng = random.Random(3)
    xs = [rng.random() for _ in range(200_000)]
    # A uniform distribution has excess kurtosis -1.2.
    assert realized_excess_kurtosis(xs) == pytest.approx(-1.2, abs=0.05)


def test_realized_skewness_of_skewed_sample():
    # Squared normals (chi-square, 1 dof) are strongly right-skewed.
    rng = random.Random(4)
    xs = [rng.gauss(0, 1) ** 2 for _ in range(100_000)]
    assert realized_skewness(xs) > 1.5


def test_moment_estimators_reject_short_series():
    with pytest.raises(ValueError):
        realized_skewness([1.0, 2.0])
    with pytest.raises(ValueError):
        realized_excess_kurtosis([1.0, 2.0, 3.0])
