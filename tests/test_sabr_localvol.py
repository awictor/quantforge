"""Tests for SABR-implied Dupire local volatility."""

import math

import pytest

from quantforge import sabr_local_vol


def test_flat_sabr_local_vol_is_alpha():
    # beta = 1, nu -> 0: the SABR smile is flat at alpha, so local vol = alpha.
    lv = sabr_local_vol(100, 100, 1.0, 0.05, alpha=0.2, beta=1.0, rho=0.0, nu=1e-6)
    assert lv == pytest.approx(0.2, abs=1e-3)


def test_positive_and_finite_across_strikes():
    for K in (80, 90, 100, 110, 120):
        lv = sabr_local_vol(100, K, 1.0, 0.05, alpha=0.25, beta=0.5, rho=-0.4, nu=0.5)
        assert lv > 0 and math.isfinite(lv)


def test_local_vol_steeper_than_implied_near_money():
    # Dupire's rule of thumb: near the money the local-vol skew is about twice
    # the implied-vol skew, so local vol at a low strike exceeds the implied.
    from quantforge import sabr_vol
    S, T, r = 100, 1.0, 0.05
    F = S * math.exp(r * T)
    K = 90
    lv = sabr_local_vol(S, K, T, r, alpha=0.25, beta=0.5, rho=-0.4, nu=0.5)
    iv = sabr_vol(F, K, T, 0.25, 0.5, -0.4, 0.5)
    assert lv > iv


def test_flat_beta_one_matches_across_strikes():
    # With nu -> 0 and beta = 1 the local vol is flat at alpha for every strike.
    for K in (85, 100, 115):
        lv = sabr_local_vol(100, K, 0.5, 0.03, alpha=0.3, beta=1.0, rho=0.0, nu=1e-6)
        assert lv == pytest.approx(0.3, abs=2e-3)
