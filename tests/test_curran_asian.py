"""Curran arithmetic-average Asian approximation."""

import math
import random

import pytest

from quantforge import curran_asian as cur
from quantforge.bsm import call_price


F, K, SIG, R, T = 100, 100, 0.25, 0.05, 1.0


def _mc(n, N, seed):
    random.seed(seed)
    disc = math.exp(-R * T)
    dt = T / n
    acc = 0.0
    for _ in range(N):
        s = F
        ssum = 0.0
        for _ in range(n):
            s *= math.exp(-0.5 * SIG * SIG * dt + SIG * math.sqrt(dt) * random.gauss(0, 1))
            ssum += s
        acc += max(ssum / n - K, 0.0)
    return disc * acc / N


def test_single_fixing_equals_black():
    assert abs(cur(F, K, SIG, R, T, 1) - call_price(F, K, T, R, SIG, b=0.0)) < 1e-3


def test_matches_monte_carlo():
    c = cur(F, K, SIG, R, T, 12)
    assert abs(c - _mc(12, 200000, 1)) < 0.05


def test_put_call_parity():
    c = cur(F, K, SIG, R, T, 12)
    p = cur(F, K, SIG, R, T, 12, is_call=False)
    assert abs((c - p) - math.exp(-R * T) * (F - K)) < 1e-6


def test_itm_above_otm():
    assert cur(F, 80, SIG, R, T, 12) > cur(F, 120, SIG, R, T, 12)


def test_validation():
    with pytest.raises(ValueError):
        cur(-1, K, SIG, R, T, 12)
    with pytest.raises(ValueError):
        cur(F, K, SIG, R, T, 0)
