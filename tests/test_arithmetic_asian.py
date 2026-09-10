"""Tests for the Turnbull-Wakeman arithmetic-average Asian closed form."""

import math

import pytest

from quantforge import (
    arithmetic_asian, geometric_asian, arithmetic_asian_mc,
    call_price, OptionType,
)


def test_arithmetic_at_least_geometric():
    # AM-GM: the arithmetic-average option is worth at least the geometric one.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.3
    a = arithmetic_asian(S, K, t, r, sigma, OptionType.CALL)
    g = geometric_asian(S, K, t, r, sigma, OptionType.CALL)
    assert a >= g


def test_arithmetic_below_vanilla():
    # Averaging cuts effective vol, so the Asian is cheaper than the vanilla.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.3
    a = arithmetic_asian(S, K, t, r, sigma, OptionType.CALL)
    assert 0 < a < call_price(S, K, t, r, sigma)


@pytest.mark.slow
@pytest.mark.parametrize("sigma", [0.15, 0.3])
@pytest.mark.parametrize("K", [90, 100, 110])
def test_matches_monte_carlo(sigma, K):
    S, t, r = 100, 1.0, 0.05
    tw = arithmetic_asian(S, K, t, r, sigma, OptionType.CALL)
    mc = arithmetic_asian_mc(S, K, t, r, sigma, OptionType.CALL,
                             n_steps=250, n_paths=100_000, control_variate=True, seed=1)
    # Turnbull-Wakeman moment-matching is accurate to a few cents vs MC.
    assert tw == pytest.approx(mc.price, abs=0.1)


def test_put_call_parity_on_average():
    # C - P = disc * (M1 - K), where M1 is the average forward.
    S, K, t, r, sigma = 100, 95, 1.0, 0.05, 0.25
    b = r
    M1 = S * (math.exp(b * t) - 1.0) / (b * t)
    c = arithmetic_asian(S, K, t, r, sigma, OptionType.CALL)
    p = arithmetic_asian(S, K, t, r, sigma, OptionType.PUT)
    assert c - p == pytest.approx(math.exp(-r * t) * (M1 - K), abs=1e-6)


def test_zero_vol_is_average_forward_intrinsic():
    S, K, t, r = 100, 90, 1.0, 0.05
    b = r
    M1 = S * (math.exp(b * t) - 1.0) / (b * t)
    v = arithmetic_asian(S, K, t, r, 0.0, OptionType.CALL)
    assert v == pytest.approx(math.exp(-r * t) * max(M1 - K, 0.0), abs=1e-9)


def test_zero_carry_finite():
    v = arithmetic_asian(100, 100, 1.0, 0.0, 0.2, OptionType.CALL, b=0.0)
    assert math.isfinite(v) and v > 0
