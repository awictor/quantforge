"""Tests for coin change / subset sum, cross-checked against brute-force references."""

import itertools
import random
from collections import Counter

import pytest

from quantforge.coin_change import min_coins, count_change, subset_sum


def _brute_min(coins, t):
    INF = float("inf")
    best = [0] + [INF] * t
    for a in range(1, t + 1):
        for c in coins:
            if c <= a:
                best[a] = min(best[a], best[a - c] + 1)
    return -1 if best[t] == INF else best[t]


def _brute_subset(vals, t):
    for r in range(len(vals) + 1):
        for combo in itertools.combinations(vals, r):
            if sum(combo) == t:
                return True
    return False


def test_fuzz_min_coins():
    rng = random.Random(181)
    for _ in range(3000):
        coins = rng.sample(range(1, 12), rng.randint(1, 5))
        t = rng.randint(0, 40)
        cnt, ms = min_coins(coins, t)
        assert cnt == _brute_min(coins, t)
        if cnt >= 0:
            assert sum(ms) == t
            assert len(ms) == cnt
            assert all(x in coins for x in ms)


def test_fuzz_count_change_ordering_independent():
    rng = random.Random(182)
    for _ in range(2000):
        coins = rng.sample(range(1, 10), rng.randint(1, 4))
        t = rng.randint(0, 30)
        # reference: brute multiset count via recursion over sorted coins
        def count(i, rem):
            if rem == 0:
                return 1
            if i >= len(coins) or rem < 0:
                return 0
            return count(i, rem - coins[i]) + count(i + 1, rem)
        assert count_change(coins, t) == count(0, t)


def test_fuzz_subset_sum():
    rng = random.Random(183)
    for _ in range(3000):
        vals = [rng.randint(0, 10) for _ in range(rng.randint(0, 10))]
        t = rng.randint(0, 30)
        ok, sub = subset_sum(vals, t)
        assert ok == _brute_subset(vals, t)
        if ok:
            assert sum(sub) == t
            assert not (Counter(sub) - Counter(vals))  # sub is a sub-multiset


def test_min_coins_explicit():
    assert min_coins([1, 5, 10, 25], 63) == (6, [1, 1, 1, 10, 25, 25])
    assert min_coins([2], 3) == (-1, [])
    assert min_coins([1, 2, 5], 0) == (0, [])


def test_min_coins_greedy_would_fail():
    # greedy picks 4+4 then can't finish 6; optimal is 3+3
    cnt, ms = min_coins([1, 3, 4], 6)
    assert cnt == 2
    assert ms == [3, 3]


def test_count_change_explicit():
    assert count_change([1, 2, 5], 5) == 4
    assert count_change([2], 3) == 0
    assert count_change([1, 2, 5], 0) == 1


def test_subset_sum_explicit():
    ok, sub = subset_sum([3, 34, 4, 12, 5, 2], 9)
    assert ok
    assert sum(sub) == 9
    assert subset_sum([1, 2, 3], 7) == (False, [])
    assert subset_sum([1, 2, 3], 0) == (True, [])


def test_subset_sum_with_zeros():
    ok, sub = subset_sum([0, 0, 5], 5)
    assert ok
    assert sum(sub) == 5


def test_negative_target_raises():
    with pytest.raises(ValueError):
        min_coins([1, 2], -1)
    with pytest.raises(ValueError):
        count_change([1, 2], -1)
    with pytest.raises(ValueError):
        subset_sum([1, 2], -1)


def test_nonpositive_coin_raises():
    with pytest.raises(ValueError):
        min_coins([0, 1], 5)
    with pytest.raises(ValueError):
        count_change([-1, 2], 5)


def test_negative_value_subset_raises():
    with pytest.raises(ValueError):
        subset_sum([1, -2, 3], 4)
