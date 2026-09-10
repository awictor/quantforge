"""Two-factor G2++ Gaussian short-rate model."""

import math
import random

import pytest

from quantforge import (
    g2pp_zero_bond,
    g2pp_bond_option,
    g2pp_caplet,
    cheyette_bond_option,
)
from quantforge.g2pp import _bond_vol


R0 = 0.03
A, B, SIGMA, ETA, RHO = 0.1, 0.3, 0.01, 0.008, -0.5


def _P0(t):
    return math.exp(-R0 * t)


def test_zero_bond_reconstitutes_curve_at_zero_state():
    P = g2pp_zero_bond(_P0(3), _P0(1), 0.0, 0.0, A, B, SIGMA, ETA, RHO, 1.0, 3.0)
    # At x = y = 0 the drift-adjustment A is tiny; bond ~ forward ratio.
    assert P == pytest.approx(_P0(3) / _P0(1), rel=1e-3)


def test_bond_option_matches_forward_measure_mc():
    e, T, K = 1.0, 3.0, 0.9
    cf = g2pp_bond_option(_P0(e), _P0(T), A, B, SIGMA, ETA, RHO, e, T, K, True)
    sig_p = _bond_vol(A, B, SIGMA, ETA, RHO, e, T)
    rng = random.Random(1)
    n = 400_000
    tot = 0.0
    fwd = _P0(T) / _P0(e)
    for _ in range(n):
        PtT = fwd * math.exp(-0.5 * sig_p * sig_p + sig_p * rng.gauss(0, 1))
        tot += max(PtT - K, 0.0)
    mc = _P0(e) * tot / n
    assert cf == pytest.approx(mc, abs=2e-4)


def test_put_call_parity():
    e, T, K = 1.0, 3.0, 0.9
    c = g2pp_bond_option(_P0(e), _P0(T), A, B, SIGMA, ETA, RHO, e, T, K, True)
    p = g2pp_bond_option(_P0(e), _P0(T), A, B, SIGMA, ETA, RHO, e, T, K, False)
    assert (c - p) == pytest.approx(_P0(T) - K * _P0(e), abs=1e-10)


def test_reduces_to_single_factor_when_second_factor_off():
    e, T, K = 1.0, 3.0, 0.9
    g2 = g2pp_bond_option(_P0(e), _P0(T), A, B, SIGMA, 1e-12, 0.0, e, T, K, True)
    chey = cheyette_bond_option(_P0(e), _P0(T), A, SIGMA, e, T, K, True)
    assert g2 == pytest.approx(chey, abs=1e-9)


def test_correlation_raises_bond_vol_and_price():
    e, T, K = 1.0, 3.0, 0.9
    vols = [_bond_vol(A, B, SIGMA, ETA, rho, e, T) for rho in (-0.9, 0.0, 0.9)]
    assert vols[0] < vols[1] < vols[2]
    calls = [g2pp_bond_option(_P0(e), _P0(T), A, B, SIGMA, ETA, rho, e, T, K, True)
             for rho in (-0.9, 0.0, 0.9)]
    assert calls[0] < calls[1] < calls[2]


def test_caplet_positive_and_bond_put_identity():
    reset, pay, strike = 1.0, 1.5, 0.03
    tau = pay - reset
    Kb = 1.0 / (1.0 + strike * tau)
    put = g2pp_bond_option(_P0(reset), _P0(pay), A, B, SIGMA, ETA, RHO,
                           reset, pay, Kb, False)
    cap = g2pp_caplet(_P0(reset), _P0(pay), A, B, SIGMA, ETA, RHO,
                      reset, pay, strike)
    assert cap > 0.0
    assert cap == pytest.approx((1.0 + strike * tau) * put, abs=1e-12)


def test_zero_vol_is_discounted_intrinsic():
    e, T, K = 1.0, 3.0, 0.9
    c = g2pp_bond_option(_P0(e), _P0(T), A, B, 1e-15, 1e-15, RHO, e, T, K, True)
    fwd = _P0(T) / _P0(e)
    assert c == pytest.approx(_P0(e) * max(fwd - K, 0.0), abs=1e-8)
