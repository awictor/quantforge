"""Tests for the vanna-volga FX smile."""

import pytest

from quantforge import VannaVolgaSmile, pillar_vols


def test_pillar_vols_from_quotes():
    p, a, c = pillar_vols(atm=0.10, rr=0.02, bf=0.005)
    assert a == pytest.approx(0.10)
    assert c == pytest.approx(0.10 + 0.005 + 0.01)   # atm + bf + rr/2
    assert p == pytest.approx(0.10 + 0.005 - 0.01)   # atm + bf - rr/2


def _smile():
    return VannaVolgaSmile(S=1.2, t=1.0, r_dom=0.02, r_for=0.01,
                           atm=0.10, rr=0.02, bf=0.005)


def test_smile_exact_at_pillars():
    sm = _smile()
    for k, v in sm.pillars():
        assert sm.vol(k) == pytest.approx(v, abs=1e-12)


def test_positive_rr_call_vol_above_put_vol():
    sm = _smile()
    pillars = sm.pillars()
    # Pillars sorted by strike: put (low K), atm, call (high K).
    assert pillars[-1][1] > pillars[0][1]


def test_butterfly_lifts_wings_above_atm():
    # With a positive butterfly the average of the two wing vols exceeds ATM.
    sm = _smile()
    (_, v_put), (_, v_atm), (_, v_call) = sm.pillars()
    assert 0.5 * (v_put + v_call) > v_atm


def test_zero_rr_zero_bf_is_flat_at_atm():
    sm = VannaVolgaSmile(S=1.2, t=1.0, r_dom=0.02, r_for=0.01,
                         atm=0.10, rr=0.0, bf=0.0)
    for k, v in sm.pillars():
        assert v == pytest.approx(0.10, abs=1e-12)
    # Interpolated vol is flat too.
    assert sm.vol(1.25) == pytest.approx(0.10, abs=1e-9)


def test_smile_interpolates_between_pillars():
    sm = _smile()
    ks = [k for k, _ in sm.pillars()]
    mid = 0.5 * (ks[0] + ks[1])
    v = sm.vol(mid)
    lo, hi = sm.pillars()[0][1], sm.pillars()[1][1]
    assert min(lo, hi) - 0.02 <= v <= max(lo, hi) + 0.02
