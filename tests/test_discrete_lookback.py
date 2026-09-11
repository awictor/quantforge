"""Discretely-monitored fixed-strike lookback (Broadie-Glasserman-Kou 1999)."""

import math

import pytest

from quantforge import discrete_fixed_strike_lookback, fixed_strike_lookback


S, K, R, SIG, T = 100.0, 100.0, 0.05, 0.2, 1.0


def _mc(n, ot, npaths=80000, seed=5):
    import random
    rng = random.Random(seed)
    disc = math.exp(-R * T)
    dt = T / n
    acc = 0.0
    for _ in range(npaths):
        logS = math.log(S)
        mx = mn = S
        for _ in range(n):
            z = rng.gauss(0.0, 1.0)
            logS += (R - 0.5 * SIG * SIG) * dt + SIG * math.sqrt(dt) * z
            s = math.exp(logS)
            mx = max(mx, s)
            mn = min(mn, s)
        acc += max(mx - K, 0.0) if ot == "call" else max(K - mn, 0.0)
    return disc * acc / npaths


@pytest.mark.slow
@pytest.mark.parametrize("ot", ["call", "put"])
def test_matches_monte_carlo(ot):
    n = 50
    cf = discrete_fixed_strike_lookback(S, K, T, R, SIG, n, ot)
    mc = _mc(n, ot)
    assert cf == pytest.approx(mc, rel=0.015)


def test_below_continuous_for_call():
    # Discrete monitoring samples the max less often -> lower call-on-max value.
    n = 50
    disc = discrete_fixed_strike_lookback(S, K, T, R, SIG, n, "call")
    cont = fixed_strike_lookback(S, K, T, R, SIG, "call")
    assert disc < cont


def test_converges_to_continuous():
    cont = fixed_strike_lookback(S, K, T, R, SIG, "call")
    coarse = discrete_fixed_strike_lookback(S, K, T, R, SIG, 12, "call")
    fine = discrete_fixed_strike_lookback(S, K, T, R, SIG, 20000, "call")
    # The BGK correction decays like 1/sqrt(n), so convergence is slow but
    # monotone: the finer grid is strictly closer to the continuous price.
    assert abs(fine - cont) < abs(coarse - cont)
    assert fine == pytest.approx(cont, rel=0.02)


def test_put_above_continuous():
    # A put on the min: discrete sampling raises the effective min -> lower min
    # captured -> the discrete put is cheaper than continuous too.
    n = 50
    disc = discrete_fixed_strike_lookback(S, K, T, R, SIG, n, "put")
    cont = fixed_strike_lookback(S, K, T, R, SIG, "put")
    assert disc < cont


def test_validation():
    with pytest.raises(ValueError):
        discrete_fixed_strike_lookback(S, K, T, R, SIG, 0, "call")  # n < 1
    with pytest.raises(ValueError):
        discrete_fixed_strike_lookback(S, -1.0, T, R, SIG, 50, "call")  # bad K
