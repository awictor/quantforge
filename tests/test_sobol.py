"""Sobol low-discrepancy sequence + Brownian-bridge QMC pricing."""

import math
import random
import statistics

import pytest

from quantforge import (
    OptionType,
    Sobol,
    brownian_bridge_path,
    sobol_european,
    sobol_asian,
    call_price,
)


def test_sobol_points_in_unit_cube():
    sob = Sobol(3)
    for _ in range(1000):
        p = sob.next()
        assert len(p) == 3
        assert all(0.0 <= x < 1.0 for x in p)


def test_sobol_first_coordinate_low_discrepancy():
    # The 1-D Sobol coordinate should cover [0,1) evenly: each quarter gets ~1/4.
    sob = Sobol(1)
    pts = [sob.next()[0] for _ in range(4096)]
    for lo in (0.0, 0.25, 0.5, 0.75):
        frac = sum(1 for x in pts if lo <= x < lo + 0.25) / len(pts)
        assert abs(frac - 0.25) < 0.02


def test_sobol_european_matches_black_scholes():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    price = sobol_european(S, K, t, r, sigma, OptionType.CALL, n_paths=8192)
    assert price == pytest.approx(call_price(S, K, t, r, sigma), abs=1e-2)


def test_sobol_converges_faster_than_pseudo_random():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    bs = call_price(S, K, t, r, sigma)
    N = 8192
    sobol_err = abs(sobol_european(S, K, t, r, sigma, n_paths=N) - bs)

    def pseudo(seed):
        rng = random.Random(seed)
        disc = math.exp(-r * t)
        tot = 0.0
        for _ in range(N):
            sT = S * math.exp((r - 0.5 * sigma * sigma) * t
                              + sigma * math.sqrt(t) * rng.gauss(0, 1))
            tot += max(sT - K, 0.0)
        return disc * tot / N

    pseudo_err = statistics.mean(abs(pseudo(s) - bs) for s in range(8))
    assert sobol_err < pseudo_err


def test_brownian_bridge_terminal_moments():
    sob = Sobol(6)
    terms = [brownian_bridge_path(sob.next(), 1.0)[-1] for _ in range(4000)]
    assert abs(statistics.mean(terms)) < 0.05
    assert statistics.pvariance(terms) == pytest.approx(1.0, abs=0.05)


def test_brownian_bridge_starts_from_zero_and_length():
    W = brownian_bridge_path([0.3, 0.6, 0.2, 0.8, 0.5, 0.1], 2.0)
    assert len(W) == 6  # W_1..W_6, W_0 = 0 implicit


def test_sobol_asian_matches_discrete_pseudo_mc():
    S, K, t, r, sigma, ns = 100, 100, 1.0, 0.05, 0.2, 6
    sob_price = sobol_asian(S, K, t, r, sigma, OptionType.CALL,
                            n_steps=ns, n_paths=16384)

    def pseudo(seed):
        rng = random.Random(seed)
        dt = t / ns
        disc = math.exp(-r * t)
        tot = 0.0
        for _ in range(40000):
            s = S
            avg = 0.0
            for _i in range(ns):
                s *= math.exp((r - 0.5 * sigma * sigma) * dt
                              + sigma * math.sqrt(dt) * rng.gauss(0, 1))
                avg += s
            tot += max(avg / ns - K, 0.0)
        return disc * tot / 40000

    ref = statistics.mean(pseudo(s) for s in range(6))
    assert sob_price == pytest.approx(ref, abs=0.05)


def test_bad_dims_raise():
    with pytest.raises(ValueError):
        Sobol(0)
    with pytest.raises(ValueError):
        sobol_asian(100, 100, 1.0, 0.05, 0.2, n_steps=99)
