"""Strike-space Greeks of the Andreasen-Huge surface (andreasen_huge_strike_greeks)."""

import pytest

from quantforge import andreasen_huge_strike_greeks


F, T = 100.0, 1.0


def _grid(lo, hi, step):
    ks = list(range(lo, hi + 1, step))
    return ks, [0.2] * len(ks)


def test_density_non_negative():
    strikes, lv = _grid(30, 200, 1)
    _ks, _dd, rnd = andreasen_huge_strike_greeks(F, strikes, T, lv)
    assert all(x >= 0.0 for x in rnd)


def test_density_integrates_to_one():
    strikes, lv = _grid(20, 300, 1)
    ks, _dd, rnd = andreasen_huge_strike_greeks(F, strikes, T, lv)
    dk = ks[1] - ks[0]
    integ = sum(x * dk for x in rnd)
    assert integ == pytest.approx(1.0, abs=0.01)


def test_density_integral_converges_with_grid_width():
    # Widening the strike range captures more tail mass -> integral -> 1.
    ks_n, lv_n = _grid(60, 140, 2)
    ks_w, lv_w = _grid(30, 200, 1)
    kn, _dn, rn = andreasen_huge_strike_greeks(F, ks_n, T, lv_n)
    kw, _dw, rw = andreasen_huge_strike_greeks(F, ks_w, T, lv_w)
    i_narrow = sum(rn) * (kn[1] - kn[0])
    i_wide = sum(rw) * (kw[1] - kw[0])
    assert i_wide > i_narrow


def test_dual_delta_monotone_in_minus_one_zero():
    strikes, lv = _grid(30, 200, 1)
    _ks, dd, _rnd = andreasen_huge_strike_greeks(F, strikes, T, lv)
    assert all(-1.0 - 1e-6 <= x <= 1e-6 for x in dd)
    assert all(dd[i] <= dd[i + 1] + 1e-9 for i in range(len(dd) - 1))


def test_bad_grid_raises():
    with pytest.raises(ValueError):
        andreasen_huge_strike_greeks(F, [90.0, 110.0], T, [0.2, 0.2])
