"""Greeks of a G2++ Bermudan swaption by CRN bumps (bermudan_swaption_g2pp_greeks)."""

import math

import pytest

from quantforge import (
    bermudan_swaption_g2pp_greeks, bermudan_swaption_g2pp,
)


def _flat(r):
    return lambda T: math.exp(-r * T)


P0 = _flat(0.03)
EX = [1.0, 2.0, 3.0]
K = 0.03
A, B, SIG, ETA, RHO = 0.1, 0.05, 0.01, 0.008, -0.7


def test_payer_d_fixed_matches_crn_fd():
    g = bermudan_swaption_g2pp_greeks(P0, EX, K, A, B, SIG, ETA, RHO,
                                      payer=True, n_paths=8000, seed=1)
    h = 1e-5
    up = bermudan_swaption_g2pp(P0, EX, K + h, A, B, SIG, ETA, RHO, True,
                                8000, seed=1)
    dn = bermudan_swaption_g2pp(P0, EX, K - h, A, B, SIG, ETA, RHO, True,
                                8000, seed=1)
    fd = (up - dn) / (2 * h)
    assert g["d_fixed"] == pytest.approx(fd, abs=1e-4)


def test_payer_d_fixed_negative():
    g = bermudan_swaption_g2pp_greeks(P0, EX, K, A, B, SIG, ETA, RHO,
                                      payer=True, n_paths=8000, seed=2)
    assert g["d_fixed"] < 0.0


def test_receiver_d_fixed_positive():
    g = bermudan_swaption_g2pp_greeks(P0, EX, K, A, B, SIG, ETA, RHO,
                                      payer=False, n_paths=8000, seed=3)
    assert g["d_fixed"] > 0.0


def test_reproducible():
    kw = dict(payer=True, n_paths=4000, seed=99)
    a = bermudan_swaption_g2pp_greeks(P0, EX, K, A, B, SIG, ETA, RHO, **kw)
    b = bermudan_swaption_g2pp_greeks(P0, EX, K, A, B, SIG, ETA, RHO, **kw)
    assert a["price"] == b["price"]
    assert a["d_fixed"] == b["d_fixed"]


def test_price_field_matches_price():
    g = bermudan_swaption_g2pp_greeks(P0, EX, K, A, B, SIG, ETA, RHO,
                                      payer=True, n_paths=6000, seed=4)
    p = bermudan_swaption_g2pp(P0, EX, K, A, B, SIG, ETA, RHO, True, 6000,
                               seed=4)
    assert g["price"] == pytest.approx(p, abs=1e-9)
