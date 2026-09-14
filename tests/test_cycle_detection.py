"""Tests for functional-graph cycle detection, cross-checked against a seen-set walk."""

import random

import pytest

from quantforge.cycle_detection import floyd_cycle, brent_cycle, cycle_elements


def _brute(f, x0):
    seen = {}
    x = x0
    i = 0
    while x not in seen:
        seen[x] = i
        x = f(x)
        i += 1
    mu = seen[x]
    return mu, i - mu


def test_fuzz_both_algorithms_vs_brute():
    rng = random.Random(301)
    for _ in range(5000):
        n = rng.randint(1, 30)
        nxt = [rng.randint(0, n - 1) for _ in range(n)]
        f = lambda x: nxt[x]
        x0 = rng.randint(0, n - 1)
        expected = _brute(f, x0)
        assert floyd_cycle(f, x0) == expected
        assert brent_cycle(f, x0) == expected


def test_fuzz_cycle_elements():
    rng = random.Random(302)
    for _ in range(2000):
        n = rng.randint(1, 30)
        nxt = [rng.randint(0, n - 1) for _ in range(n)]
        f = lambda x: nxt[x]
        x0 = rng.randint(0, n - 1)
        _, lam = _brute(f, x0)
        cyc = cycle_elements(f, x0)
        assert len(cyc) == lam
        assert len(set(cyc)) == lam
        for c in cyc:
            assert f(c) in cyc


def test_pure_cycle():
    f = lambda x: (x + 1) % 3
    assert floyd_cycle(f, 0) == (0, 3)
    assert brent_cycle(f, 0) == (0, 3)
    assert sorted(cycle_elements(f, 0)) == [0, 1, 2]


def test_rho_with_tail():
    nxt = [1, 2, 3, 4, 2]  # 0->1->2->3->4->2, tail 0,1; cycle 2,3,4
    f = lambda x: nxt[x]
    assert floyd_cycle(f, 0) == (2, 3)
    assert brent_cycle(f, 0) == (2, 3)
    assert cycle_elements(f, 0) == [2, 3, 4]


def test_self_loop():
    f = lambda x: 0
    assert floyd_cycle(f, 0) == (0, 1)
    assert brent_cycle(f, 0) == (0, 1)
    assert cycle_elements(f, 0) == [0]


def test_start_inside_cycle():
    nxt = [1, 2, 0]
    f = lambda x: nxt[x]
    assert floyd_cycle(f, 1) == (0, 3)
    assert brent_cycle(f, 1) == (0, 3)


def test_long_tail_short_cycle():
    # 0->1->...->9->9 (fixed point at 9): tail 9, cycle 1
    nxt = list(range(1, 10)) + [9]
    f = lambda x: nxt[x]
    assert floyd_cycle(f, 0) == (9, 1)
    assert brent_cycle(f, 0) == (9, 1)


def test_floyd_and_brent_agree_on_prng():
    # a small LCG modulo m is a functional graph; both must find the same period
    m = 1000
    f = lambda x: (13 * x + 7) % m
    assert floyd_cycle(f, 0) == brent_cycle(f, 0)


def test_max_iter_guard():
    # an ever-increasing map on unbounded integers never cycles
    with pytest.raises(ValueError):
        floyd_cycle(lambda x: x + 1, 0, max_iter=1000)
    with pytest.raises(ValueError):
        brent_cycle(lambda x: x + 1, 0, max_iter=1000)
