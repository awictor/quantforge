"""COS method (Fang-Oosterlee) vs Black-Scholes and the Carr-Madan pricer."""

import math

import pytest

from quantforge import (
    OptionType,
    cos_price,
    call_price,
    put_price,
    cgmy_price,
    nig_price,
)
from quantforge.cgmy import _cgmy_psi
from quantforge.nig import _nig_psi


def _gbm_psi(sigma):
    return lambda u: -0.5 * sigma * sigma * u * u


def test_cos_matches_black_scholes():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    c = cos_price(S, K, t, r, 0.0, _gbm_psi(sigma), OptionType.CALL, n_terms=256)
    assert c == pytest.approx(call_price(S, K, t, r, sigma), abs=1e-8)
    p = cos_price(S, K, t, r, 0.0, _gbm_psi(sigma), OptionType.PUT, n_terms=256)
    assert p == pytest.approx(put_price(S, K, t, r, sigma), abs=1e-8)


@pytest.mark.parametrize("K", [80, 90, 100, 110, 120])
def test_cos_matches_black_scholes_across_strikes(K):
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    c = cos_price(S, K, t, r, 0.0, _gbm_psi(sigma), OptionType.CALL, n_terms=256)
    assert c == pytest.approx(call_price(S, K, t, r, sigma), abs=1e-6)


@pytest.mark.parametrize("K", [80, 90, 100, 110, 120])
def test_cos_matches_carr_madan_cgmy(K):
    S, t, r, q = 100, 1.0, 0.03, 0.0
    C, G, M, Y = 4.0, 5.0, 10.0, 0.5
    cos = cos_price(S, K, t, r, q, lambda u: _cgmy_psi(u, C, G, M, Y),
                    OptionType.CALL, n_terms=512, L=12)
    cm = cgmy_price(S, K, t, r, C, G, M, Y, OptionType.CALL, q=q)
    assert cos == pytest.approx(cm, abs=1e-4)


@pytest.mark.parametrize("alpha,beta", [(15.0, -5.0), (20.0, -8.0), (10.0, 3.0)])
def test_cos_matches_carr_madan_nig(alpha, beta):
    S, t, r, q = 100, 1.0, 0.03, 0.0
    for K in (90, 100, 110):
        cos = cos_price(S, K, t, r, q, lambda u: _nig_psi(u, alpha, beta, 0.5),
                        OptionType.CALL, n_terms=512, L=12)
        nm = nig_price(S, K, t, r, alpha, beta, 0.5, OptionType.CALL, q=q)
        assert cos == pytest.approx(nm, abs=1e-4)


def test_cos_put_call_parity():
    S, K, t, r, q = 100, 105, 1.0, 0.04, 0.01
    psi = lambda u: _nig_psi(u, 15.0, -5.0, 0.5)
    c = cos_price(S, K, t, r, q, psi, OptionType.CALL)
    p = cos_price(S, K, t, r, q, psi, OptionType.PUT)
    rhs = S * math.exp(-q * t) - K * math.exp(-r * t)
    assert (c - p) == pytest.approx(rhs, abs=1e-4)


def test_cos_converges_in_terms():
    S, K, t, r, sigma = 100, 110, 0.5, 0.03, 0.25
    exact = call_price(S, K, t, r, sigma)
    e_lo = abs(cos_price(S, K, t, r, 0.0, _gbm_psi(sigma), n_terms=16) - exact)
    e_hi = abs(cos_price(S, K, t, r, 0.0, _gbm_psi(sigma), n_terms=256) - exact)
    assert e_hi < e_lo


def test_intrinsic_at_expiry():
    assert cos_price(120, 100, 0.0, 0.05, 0.0, _gbm_psi(0.2),
                     OptionType.CALL) == pytest.approx(20.0)
