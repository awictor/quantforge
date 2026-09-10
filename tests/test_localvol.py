"""Tests for Dupire local-volatility extraction."""

import math

import pytest

from quantforge import dupire_local_vol, local_vol_from_implied, call_price


def test_flat_implied_vol_gives_constant_local_vol():
    # A flat implied surface => local vol equals that constant everywhere.
    S, r, q = 100, 0.05, 0.02
    flat = lambda K, T: 0.25
    for K in (80, 100, 120):
        for T in (0.3, 1.0, 2.0):
            assert local_vol_from_implied(flat, S, K, T, r, q) == pytest.approx(0.25, abs=1e-4)


def test_term_structure_matches_analytic_dupire():
    # Implied total variance w(T) = 0.04 T + 0.02 T^2, strike-independent.
    # Analytic local variance = dw/dT = 0.04 + 0.04 T.
    S, r, q = 100, 0.0, 0.0
    imp = lambda K, T: math.sqrt(0.04 + 0.02 * T)
    for T in (0.5, 1.0, 2.0):
        got = local_vol_from_implied(imp, S, 100, T, r, q)
        expected = math.sqrt(0.04 + 0.04 * T)
        assert got == pytest.approx(expected, abs=1e-3)


def test_dupire_from_bsm_call_surface():
    # Build a call surface directly from BSM at constant vol; local vol = vol.
    sigma = 0.3
    call_fn = lambda K, T: call_price(100, K, T, 0.04, sigma, b=0.04)
    assert dupire_local_vol(call_fn, K=100, T=1.0, r=0.04, q=0.0) == pytest.approx(sigma, abs=1e-3)


def test_rejects_nonpositive_convexity():
    # A call price linear in K has zero convexity -> undefined local vol.
    linear = lambda K, T: max(120 - K, 0.0)
    with pytest.raises(ValueError):
        dupire_local_vol(linear, K=100, T=1.0, r=0.0)


def test_rejects_bad_inputs():
    call_fn = lambda K, T: call_price(100, K, T, 0.04, 0.2, b=0.04)
    with pytest.raises(ValueError):
        dupire_local_vol(call_fn, K=-1, T=1.0, r=0.04)
    with pytest.raises(ValueError):
        dupire_local_vol(call_fn, K=100, T=0.0, r=0.04)


def test_local_vol_positive_for_skewed_surface():
    # A mild downward skew (vol falls with strike) still yields a positive,
    # finite local vol near the money.
    S, r, q = 100, 0.03, 0.0
    skew = lambda K, T: 0.25 - 0.1 * (K / S - 1.0)
    v = local_vol_from_implied(skew, S, 100, 1.0, r, q)
    assert v > 0 and math.isfinite(v)
