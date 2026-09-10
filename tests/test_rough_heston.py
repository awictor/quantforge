"""Rough-Heston pricer via the fractional Riccati equation."""

import math

import pytest

from quantforge import (
    OptionType,
    rough_heston_price,
    rough_heston_smile,
    heston_price,
)
from quantforge.rough_heston import _rh_cf


HESTON_SETS = [
    (100, 100, 0.5, 0.03, 0.0, 0.04, 1.5, 0.04, 0.3, -0.7),
    (100, 110, 1.0, 0.05, 0.02, 0.06, 1.0, 0.05, 0.4, -0.5),
]


@pytest.mark.parametrize("S,K,t,r,q,v0,kappa,theta,nu,rho", HESTON_SETS)
def test_hurst_half_recovers_heston(S, K, t, r, q, v0, kappa, theta, nu, rho):
    # At H = 0.5 rough-Heston reduces to classical Heston with xi = kappa * nu.
    rh = rough_heston_price(S, K, t, r, v0, kappa, theta, nu, rho, H=0.5,
                            n_grid=300, q=q)
    h = heston_price(S, K, t, r, v0, kappa, theta, kappa * nu, rho,
                     OptionType.CALL, q=q)
    assert rh == pytest.approx(h, abs=5e-3)


def test_characteristic_function_is_martingale():
    # cf(-i) = E[S_T] = S e^{(r-q)t}.
    S, t, r, q = 100, 0.5, 0.03, 0.0
    cf = _rh_cf(-1j, S, t, r, q, 0.3, 1.5, 0.04, 0.3, -0.7, 0.04, 300)
    assert cf.real == pytest.approx(S * math.exp((r - q) * t), abs=1e-6)
    assert abs(cf.imag) < 1e-6


def test_put_call_parity():
    S, K, t, r, q = 100, 105, 0.5, 0.04, 0.01
    c = rough_heston_price(S, K, t, r, 0.04, 1.5, 0.04, 0.4, -0.6, H=0.2,
                           option_type=OptionType.CALL, q=q, n_grid=200)
    p = rough_heston_price(S, K, t, r, 0.04, 1.5, 0.04, 0.4, -0.6, H=0.2,
                           option_type=OptionType.PUT, q=q, n_grid=200)
    rhs = S * math.exp(-q * t) - K * math.exp(-r * t)
    assert (c - p) == pytest.approx(rhs, abs=5e-3)


def test_rough_skew_steeper_than_heston_short_maturity():
    S, r, t = 100.0, 0.0, 0.1
    strikes = [85, 92, 100, 108, 116]

    def atm_skew(H):
        sm = rough_heston_smile(S, strikes, t, r, 0.04, 1.0, 0.04, 0.4, -0.9,
                                H=H, n_grid=150)
        lm = [k for k, _ in sm]
        iv = [v for _, v in sm]
        return (iv[-1] - iv[0]) / (lm[-1] - lm[0])

    rough = atm_skew(0.1)
    classical = atm_skew(0.5)
    assert rough < 0 and classical < 0
    assert rough < classical  # steeper (more negative)


def test_intrinsic_at_expiry():
    v = rough_heston_price(100, 90, 0.0, 0.03, 0.04, 1.5, 0.04, 0.4, -0.7,
                           H=0.2, option_type=OptionType.CALL)
    assert v == pytest.approx(10.0)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        rough_heston_price(100, 100, 1.0, 0.03, 0.04, 1.5, 0.04, 0.4, -0.7, H=0.7)
    with pytest.raises(ValueError):
        rough_heston_price(100, 100, 1.0, 0.03, -0.04, 1.5, 0.04, 0.4, -0.7, H=0.2)
