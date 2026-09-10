"""Exact rainbow put closed forms (best_of_put_closed / worst_of_put_closed)."""

import math

import pytest

from quantforge import (
    best_of_put_closed,
    worst_of_put_closed,
    put_price,
)


S1, S2, K, T, R = 100.0, 95.0, 100.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


def test_put_rainbow_identity():
    # P_min + P_max = p(S1) + p(S2) at the same strike (put analogue of Stulz).
    pmin = worst_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    pmax = best_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    p1 = put_price(S1, K, T, R, SIG1)
    p2 = put_price(S2, K, T, R, SIG2)
    assert pmin + pmax == pytest.approx(p1 + p2, abs=1e-9)


def test_put_on_min_exceeds_put_on_max():
    pmin = worst_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    pmax = best_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert pmin > pmax


def test_put_on_max_below_each_vanilla_put():
    # The max is >= each asset, so a put on it is worth <= each vanilla put.
    pmax = best_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert pmax <= put_price(S1, K, T, R, SIG1) + 1e-9
    assert pmax <= put_price(S2, K, T, R, SIG2) + 1e-9


def _mc_put(kind, seed, n=600_000):
    import random
    rng = random.Random(seed)
    d1 = (R - 0.5 * SIG1 * SIG1) * T
    d2 = (R - 0.5 * SIG2 * SIG2) * T
    v1 = SIG1 * math.sqrt(T)
    v2 = SIG2 * math.sqrt(T)
    c2 = math.sqrt(1.0 - RHO * RHO)
    disc = math.exp(-R * T)
    tot = 0.0
    for _ in range(n):
        z1 = rng.gauss(0.0, 1.0)
        z2 = rng.gauss(0.0, 1.0)
        a = S1 * math.exp(d1 + v1 * z1)
        bb = S2 * math.exp(d2 + v2 * (RHO * z1 + c2 * z2))
        chosen = min(a, bb) if kind == "min" else max(a, bb)
        tot += max(K - chosen, 0.0)
    return disc * tot / n


@pytest.mark.slow
def test_put_on_min_matches_mc():
    pmin = worst_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert pmin == pytest.approx(_mc_put("min", 1), abs=0.04)


@pytest.mark.slow
def test_put_on_max_matches_mc():
    pmax = best_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert pmax == pytest.approx(_mc_put("max", 2), abs=0.04)


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        best_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, 1.5)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        worst_of_put_closed(-1, S2, K, T, R, SIG1, SIG2, RHO)
