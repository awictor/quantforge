"""Tail and shape metrics of the smile-implied risk-neutral density."""

import math

import pytest

from quantforge import (
    tail_probability,
    density_entropy,
    expected_shortfall,
)
from quantforge.mathfns import norm_cdf


S0, T, R = 100.0, 1.0, 0.05


def _flat(K):
    return 0.2


def test_upper_tail_matches_lognormal():
    F = S0 * math.exp(R * T)
    K, sig = 110.0, 0.2
    vsqrt = sig * math.sqrt(T)
    d2 = (math.log(F / K) - 0.5 * vsqrt * vsqrt) / vsqrt
    q = tail_probability(S0, T, R, _flat, 110, lower=False)
    assert q == pytest.approx(norm_cdf(d2), abs=1e-2)


def test_tails_sum_to_one():
    up = tail_probability(S0, T, R, _flat, 105, lower=False)
    dn = tail_probability(S0, T, R, _flat, 105, lower=True)
    assert (up + dn) == pytest.approx(1.0, abs=2e-2)


def test_downward_skew_fattens_left_tail():
    def down(K):
        return max(0.05, 0.2 + 0.15 * math.log(100.0 / K))
    flat_q = tail_probability(S0, T, R, _flat, 80, lower=True)
    skew_q = tail_probability(S0, T, R, down, 80, lower=True)
    assert skew_q > flat_q


def test_entropy_increases_with_vol():
    lo = density_entropy(S0, T, R, lambda K: 0.2)
    hi = density_entropy(S0, T, R, lambda K: 0.4)
    assert hi > lo


def test_expected_shortfall_in_tail():
    esd = expected_shortfall(S0, T, R, _flat, 90, lower=True)
    esu = expected_shortfall(S0, T, R, _flat, 120, lower=False)
    assert esd < 90.0
    assert esu > 120.0


def test_empty_tail_is_nan():
    # No mass below a tiny level.
    es = expected_shortfall(S0, T, R, _flat, 1.0, lower=True)
    assert math.isnan(es)
