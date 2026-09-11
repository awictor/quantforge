"""FX forward from discount curves (fxforward.fx_forward_from_curves)."""

import math

import pytest

from quantforge import fx_forward_from_curves, fx_forward, DiscountCurve


S, T = 1.10, 1.0


def _flat(r):
    return DiscountCurve.from_zero_rates([1, 2, 5], [r] * 3)


def test_matches_closed_form_on_flat_curves():
    pc, bc = _flat(0.05), _flat(0.03)
    assert fx_forward_from_curves(S, pc, bc, T) == pytest.approx(
        fx_forward(S, 0.05, 0.03, T), abs=1e-9)


def test_base_higher_rate_gives_discount():
    pc, bc = _flat(0.03), _flat(0.05)   # base rate higher
    assert fx_forward_from_curves(S, pc, bc, T) < S


def test_accepts_plain_callables():
    f = fx_forward_from_curves(
        S, lambda t: math.exp(-0.05 * t), lambda t: math.exp(-0.03 * t), T)
    assert f == pytest.approx(fx_forward(S, 0.05, 0.03, T), abs=1e-9)


def test_zero_tenor_is_spot():
    assert fx_forward_from_curves(S, _flat(0.05), _flat(0.03), 0.0) == pytest.approx(
        S, abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        fx_forward_from_curves(-1.0, _flat(0.05), _flat(0.03), T)
    with pytest.raises(ValueError):
        fx_forward_from_curves(S, _flat(0.05), _flat(0.03), -1.0)
