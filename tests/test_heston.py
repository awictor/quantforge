"""Tests for the Heston stochastic-volatility pricer.

Key checks: the xi -> 0 limit collapses to Black-Scholes (deterministic vol),
put-call parity holds, and a reference parameter set matches a value pinned
against Monte Carlo out of band.
"""

import math

import pytest

from quantforge import heston_price, call_price, put_price, OptionType
from quantforge.heston import _gauss_legendre


def test_gauss_legendre_integrates_polynomials():
    # A 64-point rule integrates x^2 on [-1,1] exactly (= 2/3).
    nodes, weights = _gauss_legendre(64)
    approx = sum(w * x * x for x, w in zip(nodes, weights))
    assert approx == pytest.approx(2.0 / 3.0, abs=1e-12)
    # And x^10 (= 2/11).
    approx10 = sum(w * x ** 10 for x, w in zip(nodes, weights))
    assert approx10 == pytest.approx(2.0 / 11.0, abs=1e-10)


@pytest.mark.parametrize("sigma", [0.15, 0.25, 0.4])
@pytest.mark.parametrize("K", [80, 100, 120])
def test_zero_volvol_reduces_to_bsm_call(sigma, K):
    # xi -> 0 with v0 = theta = sigma^2 and fast mean reversion => constant vol.
    v = sigma * sigma
    h = heston_price(100, K, 1.0, 0.05, v, 5.0, v, 1e-5, 0.0, OptionType.CALL)
    assert h == pytest.approx(call_price(100, K, 1.0, 0.05, sigma), abs=1e-3)


def test_zero_volvol_reduces_to_bsm_put():
    sigma = 0.3
    v = sigma * sigma
    h = heston_price(100, 105, 0.5, 0.03, v, 4.0, v, 1e-5, 0.0, OptionType.PUT)
    assert h == pytest.approx(put_price(100, 105, 0.5, 0.03, sigma), abs=1e-3)


def test_put_call_parity():
    kw = dict(S=100, K=95, t=1.0, r=0.04, v0=0.04, kappa=2.0, theta=0.05,
              xi=0.3, rho=-0.6, q=0.02)
    c = heston_price(**kw, option_type=OptionType.CALL)
    p = heston_price(**kw, option_type=OptionType.PUT)
    lhs = c - p
    rhs = 100 * math.exp(-0.02 * 1.0) - 95 * math.exp(-0.04 * 1.0)
    assert lhs == pytest.approx(rhs, abs=1e-6)


def test_reference_value():
    # S=K=100, t=1, r=0, v0=theta=0.04, kappa=2, xi=0.3, rho=-0.7.
    # Matches a full-truncation Euler Monte Carlo to Euler bias (~7.62).
    v = heston_price(100, 100, 1.0, 0.0, 0.04, 2.0, 0.04, 0.3, -0.7, OptionType.CALL)
    assert v == pytest.approx(7.616, abs=0.05)


def test_negative_rho_below_symmetric():
    # Negative spot/vol correlation cheapens an OTM call relative to rho=0.
    base = dict(S=100, K=115, t=1.0, r=0.0, v0=0.04, kappa=2.0, theta=0.04, xi=0.5)
    c_neg = heston_price(**base, rho=-0.7, option_type=OptionType.CALL)
    c_zero = heston_price(**base, rho=0.0, option_type=OptionType.CALL)
    assert c_neg < c_zero


def test_price_positive_and_bounded():
    c = heston_price(100, 100, 1.0, 0.05, 0.04, 2.0, 0.04, 0.4, -0.5, OptionType.CALL)
    assert 0 < c < 100


def test_zero_time_is_intrinsic():
    assert heston_price(110, 100, 0.0, 0.05, 0.04, 2, 0.04, 0.3, -0.5,
                        OptionType.CALL) == pytest.approx(10.0)
