"""Tests for Variance-Gamma pricing."""

import math

import pytest

from quantforge import variance_gamma_price as vg, call_price, put_price, OptionType


@pytest.mark.parametrize("K", [80, 100, 120])
def test_small_nu_approaches_black_scholes(K):
    # nu -> 0 and theta = 0 collapses VG to Black-Scholes.
    v = vg(100, K, 1.0, 0.05, 0.2, nu=1e-5, theta=0.0, option_type=OptionType.CALL)
    assert v == pytest.approx(call_price(100, K, 1.0, 0.05, 0.2), abs=1e-2)


def test_put_call_parity():
    kw = dict(S=100, K=95, t=1.0, r=0.05, sigma=0.2, nu=0.3, theta=-0.1)
    c = vg(**kw, option_type=OptionType.CALL)
    p = vg(**kw, option_type=OptionType.PUT)
    assert c - p == pytest.approx(100 - 95 * math.exp(-0.05), abs=1e-4)


def test_reference_matches_monte_carlo():
    # Cross-checked against a gamma-time-changed Monte Carlo (~10.45).
    v = vg(100, 100, 1.0, 0.05, 0.2, nu=0.3, theta=-0.1, option_type=OptionType.CALL)
    assert v == pytest.approx(10.45, abs=0.1)


def test_negative_theta_lifts_otm_put():
    # A negative theta creates a left skew, richening a downside put.
    sym = vg(100, 85, 1.0, 0.05, 0.2, nu=0.3, theta=0.0, option_type=OptionType.PUT)
    skew = vg(100, 85, 1.0, 0.05, 0.2, nu=0.3, theta=-0.3, option_type=OptionType.PUT)
    assert skew > sym


def test_price_positive_and_bounded():
    v = vg(100, 100, 0.5, 0.03, 0.25, nu=0.4, theta=-0.15, option_type=OptionType.CALL)
    assert 0 < v < 100


def test_zero_time_is_intrinsic():
    assert vg(110, 100, 0.0, 0.05, 0.2, nu=0.3, theta=-0.1,
              option_type=OptionType.CALL) == pytest.approx(10.0)


def test_rejects_bad_params():
    with pytest.raises(ValueError):
        vg(100, 100, 1.0, 0.05, 0.2, nu=0.0, theta=-0.1)
    with pytest.raises(ValueError):
        # Violates the martingale condition 1 - theta nu - 0.5 sigma^2 nu > 0.
        vg(100, 100, 1.0, 0.05, 0.5, nu=5.0, theta=0.5)
