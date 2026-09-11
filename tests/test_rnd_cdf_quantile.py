"""Risk-neutral CDF and quantile from an implied-vol smile."""

import math

import pytest

from quantforge import (
    risk_neutral_cdf_from_smile, risk_neutral_quantile_from_smile,
    risk_neutral_density_from_smile,
)
from quantforge.mathfns import norm_cdf


S0, T, R, SIG = 100.0, 1.0, 0.03, 0.2


def _flat(K):
    return SIG


def _bs_cdf(K):
    d2 = (math.log(S0 / K) + (R - 0.5 * SIG * SIG) * T) / (SIG * math.sqrt(T))
    return norm_cdf(-d2)


def test_flat_smile_matches_bs_cdf():
    for K in (70.0, 90.0, 100.0, 115.0, 140.0):
        assert risk_neutral_cdf_from_smile(S0, T, R, _flat, K) == pytest.approx(
            _bs_cdf(K), abs=1e-4)


def test_cdf_monotone_and_bounded():
    prev = -1.0
    for K in (50.0, 70.0, 90.0, 110.0, 130.0, 160.0, 200.0):
        F = risk_neutral_cdf_from_smile(S0, T, R, _flat, K)
        assert 0.0 <= F <= 1.0
        assert F >= prev - 1e-9
        prev = F


def test_quantile_inverts_cdf():
    for p in (0.05, 0.25, 0.5, 0.75, 0.95):
        K = risk_neutral_quantile_from_smile(S0, T, R, _flat, p)
        assert risk_neutral_cdf_from_smile(S0, T, R, _flat, K) == pytest.approx(
            p, abs=1e-6)


def test_cdf_equals_integral_of_density():
    K = 110.0
    lo, n = 1.0, 6000
    dK = (K - lo) / n
    total = 0.0
    prev = risk_neutral_density_from_smile(S0, T, R, _flat, lo)
    for i in range(1, n + 1):
        x = lo + i * dK
        cur = risk_neutral_density_from_smile(S0, T, R, _flat, x)
        total += 0.5 * (prev + cur) * dK
        prev = cur
    assert total == pytest.approx(
        risk_neutral_cdf_from_smile(S0, T, R, _flat, K), abs=2e-3)


def test_median_matches_lognormal():
    # Median of lognormal S_T = F * exp(-0.5 sigma^2 t).
    K50 = risk_neutral_quantile_from_smile(S0, T, R, _flat, 0.5)
    F = S0 * math.exp(R * T)
    assert K50 == pytest.approx(F * math.exp(-0.5 * SIG * SIG * T), rel=1e-3)


def test_quantile_rejects_out_of_range():
    with pytest.raises(ValueError):
        risk_neutral_quantile_from_smile(S0, T, R, _flat, 0.0)
    with pytest.raises(ValueError):
        risk_neutral_quantile_from_smile(S0, T, R, _flat, 1.0)


def test_skewed_smile_left_tail_heavier():
    # Equity skew (higher vol at low strikes) fattens the left tail: the 5th
    # percentile strike sits below the flat-smile 5th percentile.
    def skew(K):
        return max(0.05, SIG + 0.15 * math.log(S0 / K))  # down-strikes -> higher vol

    k_flat = risk_neutral_quantile_from_smile(S0, T, R, _flat, 0.05)
    k_skew = risk_neutral_quantile_from_smile(S0, T, R, skew, 0.05)
    assert k_skew < k_flat
