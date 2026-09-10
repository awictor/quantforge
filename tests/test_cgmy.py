"""CGMY tempered-stable Levy pricer: Carr-Madan vs Gil-Pelaez, parity, skew."""

import cmath
import math

import pytest

from quantforge import OptionType, cgmy_price, cgmy_smile
from quantforge.cgmy import _cgmy_char_logspot


# (S, K, t, r, q, (C, G, M, Y))
CASES = [
    (100, 100, 1.0, 0.03, 0.0, (4.0, 5.0, 10.0, 0.5)),
    (100, 90, 0.5, 0.05, 0.02, (2.0, 8.0, 12.0, 0.8)),
    (100, 110, 1.0, 0.04, 0.0, (1.0, 10.0, 5.0, 1.2)),
]


def _gil_pelaez_call(S, K, t, r, q, C, G, M, Y, N=2000, up=250.0):
    """Independent Gil-Pelaez inversion of the same CGMY characteristic function."""
    x0, lnK, dnu = math.log(S), math.log(K), up / N
    cf_mi = _cgmy_char_logspot(-1j, x0, t, r, q, C, G, M, Y)
    I1 = I2 = 0.0
    for k in range(1, N + 1):
        u = (k - 0.5) * dnu
        cf = _cgmy_char_logspot(u, x0, t, r, q, C, G, M, Y)
        cf1 = _cgmy_char_logspot(u - 1j, x0, t, r, q, C, G, M, Y)
        I2 += (cmath.exp(-1j * u * lnK) * cf / (1j * u)).real * dnu
        I1 += (cmath.exp(-1j * u * lnK) * cf1 / (1j * u * cf_mi)).real * dnu
    P2 = 0.5 + I2 / math.pi
    P1 = 0.5 + I1 / math.pi
    return S * math.exp(-q * t) * P1 - K * math.exp(-r * t) * P2


@pytest.mark.parametrize("S,K,t,r,q,args", CASES)
def test_carr_madan_matches_gil_pelaez(S, K, t, r, q, args):
    cm = cgmy_price(S, K, t, r, *args, OptionType.CALL, q=q, alpha=1.5)
    gp = _gil_pelaez_call(S, K, t, r, q, *args)
    assert cm == pytest.approx(gp, abs=1e-4)


@pytest.mark.parametrize("S,K,t,r,q,args", CASES)
def test_put_call_parity(S, K, t, r, q, args):
    c = cgmy_price(S, K, t, r, *args, OptionType.CALL, q=q)
    p = cgmy_price(S, K, t, r, *args, OptionType.PUT, q=q)
    rhs = S * math.exp(-q * t) - K * math.exp(-r * t)
    assert (c - p) == pytest.approx(rhs, abs=1e-4)


def test_heavier_down_tail_gives_downward_skew():
    # G < M means the down-jump tail is tempered less => richer downside vol.
    sm = cgmy_smile(100, [85, 92, 100, 108, 116], 0.25, 0.03, 3.0, 5.0, 15.0, 0.6)
    vols = [iv for _, iv in sm]
    assert vols[0] > vols[-1]


def test_symmetric_tails_give_symmetric_smile():
    sm = cgmy_smile(100, [85, 92, 100, 108, 116], 0.25, 0.03, 3.0, 10.0, 10.0, 0.6)
    vols = [iv for _, iv in sm]
    atm = vols[len(vols) // 2]
    assert vols[0] >= atm - 1e-3 and vols[-1] >= atm - 1e-3


def test_intrinsic_at_expiry():
    assert cgmy_price(100, 90, 0.0, 0.03, 4.0, 5.0, 10.0, 0.5,
                      OptionType.CALL) == pytest.approx(10.0)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        cgmy_price(100, 100, 1.0, 0.03, 4.0, 5.0, 10.0, 2.5)  # Y >= 2
    with pytest.raises(ValueError):
        cgmy_price(100, 100, 1.0, 0.03, 4.0, -5.0, 10.0, 0.5)  # G <= 0
    with pytest.raises(ValueError):
        # alpha + 1 >= M makes the Carr-Madan transform diverge.
        cgmy_price(100, 100, 1.0, 0.03, 4.0, 5.0, 2.0, 0.5, alpha=1.5)


def test_more_activity_raises_price():
    # A larger C (more jump activity) makes an ATM call more valuable.
    lo = cgmy_price(100, 100, 1.0, 0.03, 1.0, 8.0, 12.0, 0.6, OptionType.CALL)
    hi = cgmy_price(100, 100, 1.0, 0.03, 5.0, 8.0, 12.0, 0.6, OptionType.CALL)
    assert hi > lo
