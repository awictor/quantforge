"""Normal Inverse Gaussian Levy pricer via the shared Carr-Madan engine."""

import cmath
import math

import pytest

from quantforge import OptionType, nig_price, nig_smile
from quantforge.nig import _nig_psi


# (S, K, t, r, q, alpha, beta, delta)
CASES = [
    (100, 100, 1.0, 0.03, 0.0, 15.0, -5.0, 0.5),
    (100, 90, 0.5, 0.05, 0.02, 20.0, -8.0, 0.3),
    (100, 110, 1.0, 0.04, 0.0, 10.0, 3.0, 0.6),
]


def _gil_pelaez_call(S, K, t, r, q, alpha, beta, delta, N=3000, up=300.0):
    """Independent Gil-Pelaez inversion of the same NIG characteristic function."""
    x0, lnK, dnu = math.log(S), math.log(K), up / N
    omega = -_nig_psi(-1j, alpha, beta, delta)

    def cf(u):
        return cmath.exp(1j * u * (x0 + (r - q + omega) * t)
                         + t * _nig_psi(u, alpha, beta, delta))

    cfmi = cf(-1j)
    I1 = I2 = 0.0
    for k in range(1, N + 1):
        u = (k - 0.5) * dnu
        I2 += (cmath.exp(-1j * u * lnK) * cf(u) / (1j * u)).real * dnu
        I1 += (cmath.exp(-1j * u * lnK) * cf(u - 1j) / (1j * u * cfmi)).real * dnu
    return S * math.exp(-q * t) * (0.5 + I1 / math.pi) \
        - K * math.exp(-r * t) * (0.5 + I2 / math.pi)


@pytest.mark.parametrize("S,K,t,r,q,alpha,beta,delta", CASES)
def test_carr_madan_matches_gil_pelaez(S, K, t, r, q, alpha, beta, delta):
    cm = nig_price(S, K, t, r, alpha, beta, delta, OptionType.CALL, q=q)
    gp = _gil_pelaez_call(S, K, t, r, q, alpha, beta, delta)
    assert cm == pytest.approx(gp, abs=1e-4)


@pytest.mark.parametrize("S,K,t,r,q,alpha,beta,delta", CASES)
def test_put_call_parity(S, K, t, r, q, alpha, beta, delta):
    c = nig_price(S, K, t, r, alpha, beta, delta, OptionType.CALL, q=q)
    p = nig_price(S, K, t, r, alpha, beta, delta, OptionType.PUT, q=q)
    rhs = S * math.exp(-q * t) - K * math.exp(-r * t)
    assert (c - p) == pytest.approx(rhs, abs=1e-4)


def test_negative_beta_downward_skew():
    sm = nig_smile(100, [85, 92, 100, 108, 116], 0.5, 0.03, 15.0, -6.0, 0.4)
    vols = [iv for _, iv in sm]
    assert vols[0] > vols[-1]


def test_positive_beta_upward_skew():
    sm = nig_smile(100, [85, 92, 100, 108, 116], 0.5, 0.03, 15.0, 6.0, 0.4)
    vols = [iv for _, iv in sm]
    assert vols[0] < vols[-1]


def test_intrinsic_at_expiry():
    assert nig_price(100, 90, 0.0, 0.03, 15.0, -5.0, 0.5,
                     OptionType.CALL) == pytest.approx(10.0)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        nig_price(100, 100, 1.0, 0.03, -15.0, -5.0, 0.5)   # alpha <= 0
    with pytest.raises(ValueError):
        nig_price(100, 100, 1.0, 0.03, 5.0, 6.0, 0.5)      # |beta| >= alpha
    with pytest.raises(ValueError):
        # alpha^2 <= (beta + cm_alpha + 1)^2: transform diverges.
        nig_price(100, 100, 1.0, 0.03, 3.0, 1.0, 0.5, cm_alpha=1.5)


def test_larger_delta_raises_price():
    lo = nig_price(100, 100, 1.0, 0.03, 15.0, -5.0, 0.3, OptionType.CALL)
    hi = nig_price(100, 100, 1.0, 0.03, 15.0, -5.0, 0.8, OptionType.CALL)
    assert hi > lo
