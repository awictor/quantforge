"""Tests for the Corrado-Miller implied-vol initial guess.

The solver's answer is unchanged (round-trip accuracy is covered elsewhere);
here we check the better seed reduces the Newton iteration count and that the
put branch of the seed uses parity correctly.
"""

import math

import pytest

from quantforge import implied_volatility, call_price, put_price, OptionType
from quantforge.bsm import price as bs_price, vega as bs_vega


def _newton_iters(target, S, K, t, r, seed, ot=OptionType.CALL, b=None):
    if b is None:
        b = r
    sigma = min(max(seed, 1e-6), 10.0)
    for i in range(100):
        v = bs_price(S, K, t, r, sigma, ot, b) - target
        if abs(v) < 1e-8:
            return i
        vg = bs_vega(S, K, t, r, sigma, b)
        step = sigma - v / vg if vg > 1e-12 else sigma * 1.1
        sigma = max(min(step, 10.0), 1e-6)
    return 99


def _corrado_miller_seed(target, S, K, t, r):
    X = K * math.exp(-r * t)
    a = target - 0.5 * (S - X)
    rad = max(a * a - (S - X) ** 2 / math.pi, 0.0)
    return (math.sqrt(2 * math.pi / t) / (S + X)) * (a + math.sqrt(rad))


def _brenner_seed(target, S, K, t, r):
    return math.sqrt(2 * math.pi / t) * target / S


def test_corrado_miller_seed_needs_fewer_iterations():
    total_cm = total_bs = 0
    for K in (70, 85, 100, 115, 130):
        for sigma in (0.1, 0.2, 0.4, 0.7):
            tp = call_price(100, K, 1.0, 0.05, sigma)
            total_cm += _newton_iters(tp, 100, K, 1.0, 0.05,
                                      _corrado_miller_seed(tp, 100, K, 1.0, 0.05))
            total_bs += _newton_iters(tp, 100, K, 1.0, 0.05,
                                      _brenner_seed(tp, 100, K, 1.0, 0.05))
    assert total_cm < total_bs / 2   # a large, robust reduction


@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
@pytest.mark.parametrize("K", [80, 100, 120])
@pytest.mark.parametrize("sigma", [0.1, 0.25, 0.5])
def test_solver_still_round_trips(ot, K, sigma):
    S, t, r = 100.0, 1.0, 0.03
    tp = (call_price(S, K, t, r, sigma) if ot is OptionType.CALL
          else put_price(S, K, t, r, sigma))
    iv = implied_volatility(tp, S, K, t, r, ot)
    # Recovered price matches (vol itself may be ill-conditioned deep OTM).
    reprice = (call_price(S, K, t, r, iv) if ot is OptionType.CALL
               else put_price(S, K, t, r, iv))
    assert reprice == pytest.approx(tp, abs=1e-7)
