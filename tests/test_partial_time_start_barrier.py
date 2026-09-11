"""Partial-time (start) single-barrier call (Heynen-Kat 1994)."""

import math

import pytest

from quantforge import (
    partial_time_start_barrier_call, partial_time_end_barrier_call,
    barrier_option,
)
from quantforge.bsm import call_price


S, K, R, SIG, T = 100.0, 100.0, 0.05, 0.2, 1.0
H = 90.0


def test_limit_no_monitoring_approaches_vanilla():
    # t1 -> 0: barrier watched for almost no time -> vanilla call.
    pt = partial_time_start_barrier_call(S, K, H, 1e-4, T, R, SIG, "down-out")
    assert pt == pytest.approx(call_price(S, K, T, R, SIG), abs=1e-2)


def test_limit_full_monitoring_approaches_standard_barrier():
    # t1 -> T2: barrier live the whole life -> standard down-out.
    pt = partial_time_start_barrier_call(S, K, H, T - 1e-4, T, R, SIG, "down-out")
    full = barrier_option(S, K, H, T, R, SIG, "call", "down-out")
    assert pt == pytest.approx(full, abs=1e-2)


def test_between_standard_and_vanilla():
    pt = partial_time_start_barrier_call(S, K, H, 0.5, T, R, SIG, "down-out")
    full = barrier_option(S, K, H, T, R, SIG, "call", "down-out")
    van = call_price(S, K, T, R, SIG)
    assert full < pt < van


def test_monotone_decreasing_in_t1():
    vals = [partial_time_start_barrier_call(S, K, H, t1, T, R, SIG, "down-out")
            for t1 in (0.25, 0.5, 0.75)]
    assert vals[0] > vals[1] > vals[2]


def test_in_out_parity():
    ki = partial_time_start_barrier_call(S, K, H, 0.5, T, R, SIG, "down-in")
    ko = partial_time_start_barrier_call(S, K, H, 0.5, T, R, SIG, "down-out")
    assert ki + ko == pytest.approx(call_price(S, K, T, R, SIG), abs=1e-9)


@pytest.mark.slow
@pytest.mark.parametrize("t1", [0.25, 0.5, 0.75])
def test_matches_monte_carlo(t1):
    cf = partial_time_start_barrier_call(S, K, H, t1, T, R, SIG, "down-out")
    import random
    rng = random.Random(9)
    nsteps, N = 1200, 80000
    disc = math.exp(-R * T)
    dt = T / nsteps
    acc = 0.0
    for _ in range(N):
        logS = math.log(S)
        hit = False
        for i in range(nsteps):
            tcur = (i + 1) * dt
            z = rng.gauss(0.0, 1.0)
            logS += (R - 0.5 * SIG * SIG) * dt + SIG * math.sqrt(dt) * z
            if tcur <= t1 and math.exp(logS) <= H:
                hit = True
        if not hit:
            acc += max(math.exp(logS) - K, 0.0)
    mc = disc * acc / N
    assert cf == pytest.approx(mc, rel=0.02)


def test_validation():
    with pytest.raises(ValueError):
        partial_time_start_barrier_call(S, K, H, 0.0, T, R, SIG, "down-out")
    with pytest.raises(ValueError):
        partial_time_start_barrier_call(S, K, 130.0, 0.5, T, R, SIG, "up-out")
