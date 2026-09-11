"""FX forwards under covered interest parity (fxforward module)."""

import math

import pytest

from quantforge import (
    fx_forward, forward_points, fx_swap_points, implied_base_rate,
    implied_price_rate,
)


S, RP, RB, T = 1.10, 0.05, 0.03, 1.0


def test_cip_formula():
    assert fx_forward(S, RP, RB, T) == pytest.approx(
        S * math.exp((RP - RB) * T), abs=1e-12)


def test_base_premium_when_price_rate_higher():
    assert fx_forward(S, 0.05, 0.03, T) > S      # base at forward premium
    assert fx_forward(S, 0.03, 0.05, T) < S      # base at forward discount


def test_forward_points():
    assert forward_points(S, RP, RB, T) == pytest.approx(
        fx_forward(S, RP, RB, T) - S, abs=1e-12)


def test_implied_rates_round_trip():
    F = fx_forward(S, RP, RB, T)
    assert implied_base_rate(S, F, RP, T) == pytest.approx(RB, abs=1e-12)
    assert implied_price_rate(S, F, RB, T) == pytest.approx(RP, abs=1e-12)


def test_swap_points_positive_for_premium():
    assert fx_swap_points(S, RP, RB, 1.0, 2.0) > 0.0


def test_zero_tenor_is_spot():
    assert fx_forward(S, RP, RB, 0.0) == S


def test_validation():
    with pytest.raises(ValueError):
        fx_forward(-1.0, RP, RB, T)
    with pytest.raises(ValueError):
        fx_swap_points(S, RP, RB, 2.0, 1.0)      # near >= far
    with pytest.raises(ValueError):
        implied_base_rate(S, 1.12, RP, 0.0)      # t = 0
