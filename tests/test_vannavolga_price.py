"""Castagna-Mercurio vanna-volga price correction vs the market instruments."""

import math

import pytest

from quantforge import VannaVolgaSmile, OptionType
from quantforge.bsm import price as bs_price


S, T, RD, RF = 1.30, 1.0, 0.02, 0.01
ATM, RR, BF = 0.10, -0.015, 0.004


def _smile():
    return VannaVolgaSmile(S, T, RD, RF, ATM, RR, BF)


def test_reprices_market_pillars_exactly():
    vv = _smile()
    for K, sig in vv.pillars():
        mkt = bs_price(S, K, T, RD, sig, OptionType.CALL, b=RD - RF)
        assert vv.price(K, RD, RF, OptionType.CALL) == pytest.approx(mkt, abs=1e-10)


def test_price_corrected_vol_matches_pillar_vols():
    vv = _smile()
    for K, sig in vv.pillars():
        assert vv.vol_price_corrected(K, RD, RF) == pytest.approx(sig, abs=1e-8)


def test_put_call_parity_of_vv_price():
    vv = _smile()
    K = 1.35
    C = vv.price(K, RD, RF, OptionType.CALL)
    P = vv.price(K, RD, RF, OptionType.PUT)
    rhs = S * math.exp(-RF * T) - K * math.exp(-RD * T)
    assert (C - P) == pytest.approx(rhs, abs=1e-10)


def test_negative_rr_gives_put_wing_above_call_wing():
    vv = _smile()
    kp, _, kc = vv._ks
    vp = vv.vol_price_corrected(kp, RD, RF)
    vc = vv.vol_price_corrected(kc, RD, RF)
    assert vp > vc  # RR < 0 => downside vol richer


def test_intermediate_strike_vol_is_sane():
    vv = _smile()
    _, ka, kc = vv._ks
    Km = 0.5 * (ka + kc)
    vm = vv.vol_price_corrected(Km, RD, RF)
    assert 0.05 < vm < 0.20


def test_flat_market_gives_flat_smile():
    # No skew, no convexity -> the correction vanishes and the smile is flat.
    vv = VannaVolgaSmile(S, T, RD, RF, atm=0.10, rr=0.0, bf=0.0)
    for K in (1.20, 1.30, 1.42):
        assert vv.vol_price_corrected(K, RD, RF) == pytest.approx(0.10, abs=1e-6)
