"""Andersen QE Heston Monte Carlo vs the Fourier char-function price."""

import math

import pytest

from quantforge import OptionType, heston_price, heston_qe_mc
from quantforge.heston_mc import _norm_ppf


# (S, K, t, r, v0, kappa, theta, xi, rho, q)
CASES = [
    (100, 100, 1.0, 0.0, 0.04, 1.0, 0.04, 0.5, -0.7, 0.0),
    (100, 100, 5.0, 0.05, 0.09, 1.0, 0.09, 1.0, -0.3, 0.02),  # long, Feller-violating
    (100, 90, 1.0, 0.03, 0.04, 2.0, 0.04, 0.6, -0.6, 0.0),
    (100, 110, 2.0, 0.04, 0.05, 1.2, 0.05, 0.7, -0.5, 0.01),
]


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,v0,kappa,theta,xi,rho,q", CASES)
def test_qe_mc_matches_char_function(S, K, t, r, v0, kappa, theta, xi, rho, q):
    cf = heston_price(S, K, t, r, v0, kappa, theta, xi, rho,
                      OptionType.CALL, q=q)
    mc = heston_qe_mc(S, K, t, r, v0, kappa, theta, xi, rho,
                      OptionType.CALL, q=q, n_steps=100, n_paths=120_000, seed=5)
    # Within ~3.5 standard errors of the exact Fourier price.
    assert abs(mc.price - cf) < 3.5 * mc.std_error + 0.02


@pytest.mark.slow
def test_qe_mc_put_matches():
    S, K, t, r, q = 100, 110, 1.0, 0.04, 0.01
    cf = heston_price(S, K, t, r, 0.04, 1.2, 0.05, 0.6, -0.6,
                      OptionType.PUT, q=q)
    mc = heston_qe_mc(S, K, t, r, 0.04, 1.2, 0.05, 0.6, -0.6,
                      OptionType.PUT, q=q, n_steps=100, n_paths=120_000, seed=2)
    assert abs(mc.price - cf) < 3.5 * mc.std_error + 0.02


@pytest.mark.slow
def test_qe_mc_is_martingale():
    # Discounted spot is a martingale: E[S_T] e^{-rt} = S e^{-qt}. Price a call
    # struck at ~0 to recover the forward.
    S, r, q, t = 100, 0.05, 0.02, 2.0
    mc = heston_qe_mc(S, 1e-6, t, r, 0.04, 1.5, 0.04, 0.6, -0.7,
                      OptionType.CALL, q=q, n_steps=100, n_paths=200_000, seed=9)
    fwd_mc = mc.price * math.exp(r * t)
    fwd_exact = S * math.exp((r - q) * t)
    assert abs(fwd_mc - fwd_exact) / fwd_exact < 1e-3


def test_norm_ppf_reference_values():
    assert _norm_ppf(0.5) == pytest.approx(0.0, abs=1e-9)
    assert _norm_ppf(0.975) == pytest.approx(1.959964, abs=1e-4)
    assert _norm_ppf(0.025) == pytest.approx(-1.959964, abs=1e-4)


def test_variance_stays_nonnegative_smoke():
    # A tiny, cheap run: QE must never blow up even with a huge vol-of-vol.
    mc = heston_qe_mc(100, 100, 1.0, 0.0, 0.04, 0.3, 0.04, 2.0, -0.9,
                      n_steps=50, n_paths=2_000, seed=1)
    assert mc.price > 0.0 and math.isfinite(mc.price)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        heston_qe_mc(100, 100, 1.0, 0.0, -0.04, 1.0, 0.04, 0.5, -0.7, n_paths=100)
    with pytest.raises(ValueError):
        heston_qe_mc(-1, 100, 1.0, 0.0, 0.04, 1.0, 0.04, 0.5, -0.7, n_paths=100)
