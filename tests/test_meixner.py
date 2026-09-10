"""Meixner Levy pricer: Carr-Madan vs COS, parity, skew."""

import math

import pytest

from quantforge import OptionType, meixner_price, meixner_smile
from quantforge.meixner import _meixner_psi
from quantforge.carrmadan import cos_price


# (a, b, d)
PARAM_SETS = [(0.3, -0.5, 0.5), (0.2, -0.3, 1.0), (0.25, 0.4, 0.8)]


@pytest.mark.parametrize("a,b,d", PARAM_SETS)
@pytest.mark.parametrize("K", [90, 100, 110])
def test_carr_madan_matches_cos(a, b, d, K):
    S, t, r, q = 100, 1.0, 0.03, 0.0
    cm = meixner_price(S, K, t, r, a, b, d, OptionType.CALL, q=q)
    cos = cos_price(S, K, t, r, q, lambda u: _meixner_psi(u, a, b, d),
                    OptionType.CALL, n_terms=512, L=12)
    assert cm == pytest.approx(cos, abs=1e-4)


@pytest.mark.parametrize("a,b,d", PARAM_SETS)
def test_put_call_parity(a, b, d):
    S, K, t, r, q = 100, 105, 1.0, 0.04, 0.01
    c = meixner_price(S, K, t, r, a, b, d, OptionType.CALL, q=q)
    p = meixner_price(S, K, t, r, a, b, d, OptionType.PUT, q=q)
    rhs = S * math.exp(-q * t) - K * math.exp(-r * t)
    assert (c - p) == pytest.approx(rhs, abs=1e-4)


def test_negative_b_downward_skew():
    sm = meixner_smile(100, [85, 92, 100, 108, 116], 0.5, 0.03, 0.3, -0.6, 0.5)
    vols = [iv for _, iv in sm]
    assert vols[0] > vols[-1]


def test_positive_b_upward_skew():
    sm = meixner_smile(100, [85, 92, 100, 108, 116], 0.5, 0.03, 0.3, 0.6, 0.5)
    vols = [iv for _, iv in sm]
    assert vols[0] < vols[-1]


def test_intrinsic_at_expiry():
    assert meixner_price(100, 90, 0.0, 0.03, 0.3, -0.5, 0.5,
                         OptionType.CALL) == pytest.approx(10.0)


def test_more_activity_raises_price():
    lo = meixner_price(100, 100, 1.0, 0.03, 0.3, -0.3, 0.3, OptionType.CALL)
    hi = meixner_price(100, 100, 1.0, 0.03, 0.3, -0.3, 1.0, OptionType.CALL)
    assert hi > lo


def test_bad_params_raise():
    with pytest.raises(ValueError):
        meixner_price(100, 100, 1.0, 0.03, -0.3, -0.5, 0.5)   # a <= 0
    with pytest.raises(ValueError):
        meixner_price(100, 100, 1.0, 0.03, 0.3, 4.0, 0.5)     # b out of (-pi, pi)
    with pytest.raises(ValueError):
        # a (cm_alpha + 1) + b too close to pi -> transform diverges.
        meixner_price(100, 100, 1.0, 0.03, 1.4, 0.5, 0.5, cm_alpha=1.5)
