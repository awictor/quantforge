"""Powered option: payoff max(S_T-K,0)**p / max(K-S_T,0)**p (exotics.powered_option)."""

import math

import pytest

from quantforge import powered_option, powered_option_greeks
from quantforge.bsm import price, OptionType


def _mc(S, K, t, r, sig, p, ot, n=400000, seed=0):
    import random
    rng = random.Random(seed)
    acc = 0.0
    for _ in range(n):
        z = rng.gauss(0.0, 1.0)
        for zz in (z, -z):  # antithetic
            ST = S * math.exp((r - 0.5 * sig * sig) * t + sig * math.sqrt(t) * zz)
            pay = max(ST - K, 0.0) if ot == "call" else max(K - ST, 0.0)
            acc += pay ** p
    return math.exp(-r * t) * acc / (2 * n)


def test_power_one_recovers_vanilla():
    S, K, t, r, sig = 100.0, 105.0, 1.0, 0.05, 0.2
    for ot in (OptionType.CALL, OptionType.PUT):
        assert powered_option(S, K, t, r, sig, 1, ot) == pytest.approx(
            price(S, K, t, r, sig, ot), abs=1e-10)


@pytest.mark.parametrize("ot", ["call", "put"])
def test_matches_monte_carlo_p2(ot):
    S, K, t, r, sig, p = 100.0, 100.0, 1.0, 0.05, 0.2, 2
    cf = powered_option(S, K, t, r, sig, p, ot)
    mc = _mc(S, K, t, r, sig, p, ot, n=400000, seed=1)
    assert cf == pytest.approx(mc, rel=0.02)


def test_matches_monte_carlo_p3_call():
    S, K, t, r, sig, p = 100.0, 95.0, 0.75, 0.03, 0.25, 3
    cf = powered_option(S, K, t, r, sig, p, "call")
    mc = _mc(S, K, t, r, sig, p, "call", n=400000, seed=2)
    assert cf == pytest.approx(mc, rel=0.03)


def test_zero_vol_is_discounted_intrinsic_power():
    S, K, t, r = 100.0, 90.0, 1.0, 0.05
    fwd = S * math.exp(r * t)
    expect = math.exp(-r * t) * max(fwd - K, 0.0) ** 2
    assert powered_option(S, K, t, r, 0.0, 2, "call") == pytest.approx(expect, abs=1e-9)


def test_greeks_delta_gamma_signs_call():
    g = powered_option_greeks(100.0, 100.0, 1.0, 0.05, 0.2, 2, "call")
    assert g["delta"] > 0.0   # call rises with spot
    assert g["gamma"] > 0.0   # convex payoff
    assert g["vega"] > 0.0    # more vol -> more value
    assert g["price"] == pytest.approx(
        powered_option(100.0, 100.0, 1.0, 0.05, 0.2, 2, "call"), abs=1e-9)


def test_greeks_reduce_to_bsm_at_power_one():
    from quantforge.bsm import delta as bs_delta, gamma as bs_gamma, vega as bs_vega
    S, K, t, r, sig = 100.0, 100.0, 1.0, 0.05, 0.2
    g = powered_option_greeks(S, K, t, r, sig, 1, "call")
    assert g["delta"] == pytest.approx(bs_delta(S, K, t, r, sig, "call"), abs=1e-5)
    assert g["gamma"] == pytest.approx(bs_gamma(S, K, t, r, sig), abs=1e-4)
    assert g["vega"] == pytest.approx(bs_vega(S, K, t, r, sig), rel=1e-4)


def test_requires_positive_integer_power():
    with pytest.raises(ValueError):
        powered_option(100.0, 100.0, 1.0, 0.05, 0.2, 0, "call")
    with pytest.raises(ValueError):
        powered_option(100.0, 100.0, 1.0, 0.05, 0.2, 2.5, "call")
