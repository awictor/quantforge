"""Variance-Gamma on the shared Carr-Madan engine: COS cross-check + smile."""

import math

import pytest

from quantforge import (
    OptionType,
    variance_gamma_price as vg,
    variance_gamma_smile,
    call_price,
)
from quantforge.variancegamma import _vg_psi
from quantforge.carrmadan import cos_price


PARAM_SETS = [(0.2, 0.3, -0.2), (0.15, 0.5, -0.3), (0.25, 0.2, 0.1)]


@pytest.mark.parametrize("sigma,nu,theta", PARAM_SETS)
@pytest.mark.parametrize("K", [90, 100, 110])
def test_carr_madan_matches_cos(sigma, nu, theta, K):
    S, t, r, q = 100, 1.0, 0.03, 0.0
    cm = vg(S, K, t, r, sigma, nu, theta, OptionType.CALL, q=q)
    cos = cos_price(S, K, t, r, q, lambda u: _vg_psi(u, sigma, nu, theta),
                    OptionType.CALL, n_terms=512, L=12)
    assert cm == pytest.approx(cos, abs=1e-4)


def test_small_nu_approaches_black_scholes():
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    for K in (90, 100, 110):
        v = vg(S, K, t, r, sigma, 1e-4, 0.0, OptionType.CALL)
        assert v == pytest.approx(call_price(S, K, t, r, sigma), abs=1e-2)


def test_put_call_parity():
    S, K, t, r, q = 100, 105, 1.0, 0.04, 0.01
    c = vg(S, K, t, r, 0.2, 0.3, -0.2, OptionType.CALL, q=q)
    p = vg(S, K, t, r, 0.2, 0.3, -0.2, OptionType.PUT, q=q)
    rhs = S * math.exp(-q * t) - K * math.exp(-r * t)
    assert (c - p) == pytest.approx(rhs, abs=1e-4)


def test_negative_theta_downward_skew():
    sm = variance_gamma_smile(100, [85, 92, 100, 108, 116], 0.5, 0.03,
                              0.2, 0.3, -0.3)
    vols = [iv for _, iv in sm]
    assert vols[0] > vols[-1]


def test_rejects_bad_params():
    with pytest.raises(ValueError):
        vg(100, 100, 1.0, 0.03, 0.2, -0.1, -0.2)   # nu <= 0
    with pytest.raises(ValueError):
        vg(100, 100, 1.0, 0.03, 0.2, 5.0, 0.3)     # violates martingale cond.
