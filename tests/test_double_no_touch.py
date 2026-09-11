"""Double-no-touch / double-one-touch binaries (double_no_touch, double_one_touch)."""

import math

import pytest

from quantforge import (
    double_no_touch, double_one_touch, no_touch,
)


S, T, R = 100.0, 0.5, 0.03


def test_upper_barrier_to_infinity_matches_lower_no_touch():
    # As U -> infinity only the lower barrier can knock: DNT -> no_touch(L).
    for L, sig in ((90.0, 0.2), (80.0, 0.25), (95.0, 0.15)):
        dnt = double_no_touch(S, L, 1e7, T, R, sig)
        assert dnt == pytest.approx(no_touch(S, L, T, R, sig), rel=1e-5)


def test_lower_barrier_to_zero_matches_upper_no_touch():
    for U, sig in ((115.0, 0.2), (130.0, 0.25), (108.0, 0.15)):
        dnt = double_no_touch(S, 1e-3, U, T, R, sig)
        assert dnt == pytest.approx(no_touch(S, U, T, R, sig), rel=1e-6)


def test_dnt_plus_dot_equals_discounted_cash():
    disc = math.exp(-R * T)
    for L, U, sig in ((90.0, 115.0, 0.2), (85.0, 130.0, 0.25)):
        d = double_no_touch(S, L, U, T, R, sig)
        o = double_one_touch(S, L, U, T, R, sig)
        assert d + o == pytest.approx(disc, abs=1e-9)


def test_narrower_band_lower_survival():
    wide = double_no_touch(S, 85.0, 120.0, T, R, 0.2)
    narrow = double_no_touch(S, 95.0, 108.0, T, R, 0.2)
    assert narrow < wide


def test_survival_between_zero_and_discounted_cash():
    disc = math.exp(-R * T)
    v = double_no_touch(S, 90.0, 115.0, T, R, 0.2)
    assert 0.0 <= v <= disc


def test_higher_vol_lowers_survival():
    lo = double_no_touch(S, 88.0, 114.0, T, R, 0.15)
    hi = double_no_touch(S, 88.0, 114.0, T, R, 0.35)
    assert hi < lo


@pytest.mark.slow
def test_matches_fine_step_monte_carlo():
    import random

    L, U, sig = 88.0, 116.0, 0.2
    closed = double_no_touch(S, L, U, T, R, sig)
    # Continuously-monitored survival: refine the grid so discrete-monitoring
    # bias (which overstates survival) is small, then allow a modest tolerance.
    rng = random.Random(12345)
    steps, paths = 4000, 40000
    dt = T / steps
    drift = (R - 0.5 * sig * sig) * dt
    vol = sig * math.sqrt(dt)
    surv = 0
    for _ in range(paths):
        x = S
        alive = True
        for _ in range(steps):
            x *= math.exp(drift + vol * rng.gauss(0.0, 1.0))
            if x <= L or x >= U:
                alive = False
                break
        if alive:
            surv += 1
    mc = math.exp(-R * T) * surv / paths
    assert mc == pytest.approx(closed, abs=0.02)


def test_requires_ordering():
    with pytest.raises(ValueError):
        double_no_touch(S, 110.0, 90.0, T, R, 0.2)   # L > U
    with pytest.raises(ValueError):
        double_no_touch(S, 100.0, 120.0, T, R, 0.2)  # S == L
