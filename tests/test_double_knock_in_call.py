"""Double-barrier knock-in call via in-out parity (double_knock_in_call)."""

import pytest

from quantforge import double_knock_in_call, double_knock_out_call
from quantforge.bsm import call_price


S, K, T, R, SIG = 100.0, 95.0, 0.5, 0.05, 0.25


@pytest.mark.parametrize("L,U", [(90.0, 115.0), (80.0, 130.0), (85.0, 120.0)])
def test_in_out_parity(L, U):
    dki = double_knock_in_call(S, K, L, U, T, R, SIG)
    dko = double_knock_out_call(S, K, L, U, T, R, SIG)
    assert dki + dko == pytest.approx(call_price(S, K, T, R, SIG), abs=1e-9)


def test_wide_barriers_near_zero():
    # Barriers far away are almost never touched: the knock-in rarely activates.
    v = double_knock_in_call(S, K, 40.0, 400.0, T, R, SIG)
    assert v == pytest.approx(0.0, abs=1e-6)


def test_tighter_band_raises_value():
    # A tighter corridor knocks in more often -> higher knock-in value.
    wide = double_knock_in_call(S, K, 80.0, 130.0, T, R, SIG)
    narrow = double_knock_in_call(S, K, 90.0, 112.0, T, R, SIG)
    assert narrow > wide


def test_nonnegative_and_below_vanilla():
    v = double_knock_in_call(S, K, 88.0, 116.0, T, R, SIG)
    assert 0.0 <= v <= call_price(S, K, T, R, SIG)


def test_requires_ordering():
    with pytest.raises(ValueError):
        double_knock_in_call(S, K, 120.0, 90.0, T, R, SIG)
