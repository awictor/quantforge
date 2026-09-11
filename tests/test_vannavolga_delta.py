"""Delta-space queries on a vanna-volga FX smile (vol_at_delta, risk_reversal, butterfly)."""

import pytest

from quantforge import VannaVolgaSmile


S, T, RD, RF = 1.30, 0.5, 0.02, 0.01
ATM, RR, BF = 0.10, -0.015, 0.004


def _smile(call_delta=0.25):
    return VannaVolgaSmile(S, T, RD, RF, ATM, RR, BF, call_delta=call_delta)


def test_vol_at_pillar_delta_matches_pillar_vols():
    sm = _smile()
    s_p, s_a, s_c = ATM + BF - 0.5 * RR, ATM, ATM + BF + 0.5 * RR
    assert sm.vol_at_delta(0.25, True) == pytest.approx(s_c, abs=1e-6)
    assert sm.vol_at_delta(0.25, False) == pytest.approx(s_p, abs=1e-6)


def test_risk_reversal_recovers_input():
    assert _smile().risk_reversal(0.25) == pytest.approx(RR, abs=1e-6)


def test_butterfly_recovers_input():
    assert _smile().butterfly(0.25) == pytest.approx(BF, abs=1e-6)


def test_atm_delta_near_atm_vol():
    # A ~50-delta option sits near the ATM pillar vol.
    sm = _smile()
    v = sm.vol_at_delta(0.5, True)
    assert v == pytest.approx(ATM, abs=5e-3)


def test_negative_rr_puts_richer_than_calls():
    # RR < 0 => 25d put vol above 25d call vol.
    sm = _smile()
    assert sm.vol_at_delta(0.25, False) > sm.vol_at_delta(0.25, True)


def test_delta_out_of_range_raises():
    sm = _smile()
    with pytest.raises(ValueError):
        sm.vol_at_delta(0.0, True)
    with pytest.raises(ValueError):
        sm.vol_at_delta(1.0, False)


def test_ten_delta_wider_than_twentyfive():
    # A convex (positive-butterfly) smile is higher in the 10d wings than the
    # 25d wings: the butterfly measured at 10 delta exceeds the 25-delta one.
    sm = _smile()
    assert sm.butterfly(0.10) > sm.butterfly(0.25)
