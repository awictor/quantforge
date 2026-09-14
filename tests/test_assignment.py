"""Combinatorial optimization: Hungarian assignment and 0/1 knapsack."""

import itertools
import random

import pytest

from quantforge import hungarian, knapsack_01


def test_hungarian_matches_brute_force():
    def brute(cost):
        n = len(cost)
        return min(sum(cost[i][perm[i]] for i in range(n))
                   for perm in itertools.permutations(range(n)))

    rng = random.Random(1)
    for _ in range(300):
        n = rng.randint(1, 5)
        cost = [[rng.randint(0, 20) for _ in range(n)] for _ in range(n)]
        a, tot = hungarian(cost)
        assert sorted(a) == list(range(n))       # valid permutation
        assert tot == brute(cost)


def test_hungarian_known():
    a, tot = hungarian([[4, 1, 3], [2, 0, 5], [3, 2, 2]])
    assert tot == 5
    assert sorted(a) == [0, 1, 2]


def test_knapsack_matches_brute_force():
    def brute(w, v, cap):
        n = len(w)
        best = 0
        for r in range(n + 1):
            for combo in itertools.combinations(range(n), r):
                if sum(w[i] for i in combo) <= cap:
                    best = max(best, sum(v[i] for i in combo))
        return best

    rng = random.Random(2)
    for _ in range(500):
        n = rng.randint(0, 8)
        w = [rng.randint(1, 10) for _ in range(n)]
        v = [rng.randint(1, 20) for _ in range(n)]
        cap = rng.randint(0, 25)
        val, chosen = knapsack_01(w, v, cap)
        assert val == brute(w, v, cap)
        assert sum(w[i] for i in chosen) <= cap
        assert sum(v[i] for i in chosen) == val   # trace is consistent


def test_knapsack_known():
    val, chosen = knapsack_01([2, 3, 4, 5], [3, 4, 5, 6], 5)
    assert val == 7
    assert sorted(chosen) == [0, 1]


def test_validation():
    with pytest.raises(ValueError):
        hungarian([[1, 2]])                       # non-square
    with pytest.raises(ValueError):
        knapsack_01([1, 2], [1], 5)               # length mismatch
    with pytest.raises(ValueError):
        knapsack_01([1], [1], -1)                 # negative capacity
