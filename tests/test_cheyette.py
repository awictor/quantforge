"""Single-factor Cheyette (quasi-Gaussian) short-rate model."""

import math
import random

import pytest

from quantforge import (
    cheyette_G,
    cheyette_y,
    cheyette_zero_bond,
    cheyette_bond_option,
    cheyette_caplet,
)


R0, KAPPA, SIGMA = 0.03, 0.1, 0.01


def _P0(t):
    return math.exp(-R0 * t)   # flat-curve discount factor


def test_G_function_limits():
    # G(tau) -> tau as kappa -> 0, and is increasing in tau.
    assert cheyette_G(1e-13, 2.0) == pytest.approx(2.0, abs=1e-6)
    assert cheyette_G(0.1, 3.0) > cheyette_G(0.1, 1.0)


def test_zero_bond_reconstitutes_curve_at_zero_state():
    # With x = y = 0 the bond is just the forward discount ratio.
    P = cheyette_zero_bond(_P0(3), _P0(1), 0.0, 0.0, KAPPA, 1.0, 3.0)
    assert P == pytest.approx(_P0(3) / _P0(1))


def test_bond_option_matches_forward_measure_mc():
    expiry, maturity, strike = 1.0, 3.0, 0.9
    cf = cheyette_bond_option(_P0(expiry), _P0(maturity), KAPPA, SIGMA,
                              expiry, maturity, strike, is_call=True)
    var = SIGMA * SIGMA * (1 - math.exp(-2 * KAPPA * expiry)) / (2 * KAPPA)
    G = cheyette_G(KAPPA, maturity - expiry)
    sig_p = G * math.sqrt(var)
    rng = random.Random(1)
    n = 400_000
    tot = 0.0
    fwd = _P0(maturity) / _P0(expiry)
    for _ in range(n):
        PtT = fwd * math.exp(-0.5 * sig_p * sig_p + sig_p * rng.gauss(0, 1))
        tot += max(PtT - strike, 0.0)
    mc = _P0(expiry) * tot / n
    assert cf == pytest.approx(mc, abs=2e-4)


def test_bond_option_put_call_parity():
    e, T, K = 1.0, 3.0, 0.9
    c = cheyette_bond_option(_P0(e), _P0(T), KAPPA, SIGMA, e, T, K, is_call=True)
    p = cheyette_bond_option(_P0(e), _P0(T), KAPPA, SIGMA, e, T, K, is_call=False)
    assert (c - p) == pytest.approx(_P0(T) - K * _P0(e), abs=1e-10)


def test_caplet_equals_bond_put_identity():
    reset, pay, strike = 1.0, 1.5, 0.03
    tau = pay - reset
    Kb = 1.0 / (1.0 + strike * tau)
    put = cheyette_bond_option(_P0(reset), _P0(pay), KAPPA, SIGMA,
                               reset, pay, Kb, is_call=False)
    cap = cheyette_caplet(_P0(reset), _P0(pay), KAPPA, SIGMA, reset, pay, strike)
    assert cap == pytest.approx((1.0 + strike * tau) * put, abs=1e-12)


def test_zero_vol_bond_option_is_discounted_intrinsic():
    e, T, K = 1.0, 3.0, 0.9
    c = cheyette_bond_option(_P0(e), _P0(T), KAPPA, 1e-15, e, T, K, is_call=True)
    fwd = _P0(T) / _P0(e)
    assert c == pytest.approx(_P0(e) * max(fwd - K, 0.0), abs=1e-8)


def test_caplet_positive():
    assert cheyette_caplet(_P0(1.0), _P0(1.5), KAPPA, SIGMA, 1.0, 1.5, 0.03) > 0.0
