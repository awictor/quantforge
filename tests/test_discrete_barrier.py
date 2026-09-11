"""Discretely-monitored single-barrier options (Broadie-Glasserman-Kou 1999)."""

import math

import pytest

from quantforge import discrete_barrier_option, barrier_option
from quantforge.exotics import Barrier


S, K, R, SIG, T = 100.0, 100.0, 0.05, 0.2, 1.0


def _mc(n, H, kind, npaths=80000, seed=5):
    import random
    rng = random.Random(seed)
    disc = math.exp(-R * T)
    dt = T / n
    down = "down" in kind
    out = "out" in kind
    acc = 0.0
    for _ in range(npaths):
        logS = math.log(S)
        hit = False
        for _ in range(n):
            z = rng.gauss(0.0, 1.0)
            logS += (R - 0.5 * SIG * SIG) * dt + SIG * math.sqrt(dt) * z
            s = math.exp(logS)
            if (down and s <= H) or (not down and s >= H):
                hit = True
        ST = math.exp(logS)
        payoff = max(ST - K, 0.0)
        alive = (not hit) if out else hit
        if alive:
            acc += payoff
    return disc * acc / npaths


@pytest.mark.slow
@pytest.mark.parametrize("H,kind", [
    (90.0, "down-out"),
    (130.0, "up-out"),
    (90.0, "down-in"),
    (130.0, "up-in"),
])
def test_matches_monte_carlo(H, kind):
    n = 50
    cf = discrete_barrier_option(S, K, H, T, R, SIG, n, "call", kind)
    mc = _mc(n, H, kind)
    assert cf == pytest.approx(mc, rel=0.03)


def test_knockout_above_continuous():
    # Discrete monitoring breaches less often -> knock-out worth more.
    n = 50
    disc = discrete_barrier_option(S, K, 90.0, T, R, SIG, n, "call", "down-out")
    cont = barrier_option(S, K, 90.0, T, R, SIG, "call", "down-out")
    assert disc > cont


def test_knockin_below_continuous():
    n = 50
    disc = discrete_barrier_option(S, K, 90.0, T, R, SIG, n, "call", "down-in")
    cont = barrier_option(S, K, 90.0, T, R, SIG, "call", "down-in")
    assert disc < cont


def test_in_out_parity_at_shifted_barrier():
    # KI + KO with the same monitoring must equal the vanilla priced at the
    # shifted barrier (both legs use the identical H_eff), i.e. the vanilla.
    from quantforge.bsm import call_price
    n = 50
    ki = discrete_barrier_option(S, K, 90.0, T, R, SIG, n, "call", "down-in")
    ko = discrete_barrier_option(S, K, 90.0, T, R, SIG, n, "call", "down-out")
    assert ki + ko == pytest.approx(call_price(S, K, T, R, SIG), abs=1e-9)


def test_converges_to_continuous():
    cont = barrier_option(S, K, 90.0, T, R, SIG, "call", "down-out")
    coarse = discrete_barrier_option(S, K, 90.0, T, R, SIG, 12, "call", "down-out")
    fine = discrete_barrier_option(S, K, 90.0, T, R, SIG, 20000, "call", "down-out")
    assert abs(fine - cont) < abs(coarse - cont)
    assert fine == pytest.approx(cont, rel=0.02)


def test_validation():
    with pytest.raises(ValueError):
        discrete_barrier_option(S, K, 90.0, T, R, SIG, 0, "call", "down-out")
    with pytest.raises(ValueError):
        discrete_barrier_option(S, K, -1.0, T, R, SIG, 50, "call", "down-out")
