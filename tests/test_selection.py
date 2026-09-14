"""Selection: k-th order statistic, median, top-k without a full sort."""

import random
import statistics

import pytest

from quantforge import kth_smallest, median, top_k


def test_kth_smallest_with_duplicates():
    rng = random.Random(1)
    for _ in range(500):
        n = rng.randint(1, 60)
        a = [rng.randint(0, 10) for _ in range(n)]     # many duplicates
        s = sorted(a)
        for k in range(n):
            assert kth_smallest(a, k) == s[k]


def test_kth_smallest_distinct():
    rng = random.Random(2)
    for _ in range(300):
        a = random.sample(range(1000), rng.randint(1, 50))
        s = sorted(a)
        for k in range(len(a)):
            assert kth_smallest(a, k) == s[k]


def test_median_matches_statistics():
    rng = random.Random(3)
    for _ in range(500):
        a = [rng.uniform(-100, 100) for _ in range(rng.randint(1, 40))]
        assert abs(median(a) - statistics.median(a)) < 1e-9


def test_top_k():
    rng = random.Random(4)
    for _ in range(300):
        n = rng.randint(1, 50)
        a = [rng.randint(0, 100) for _ in range(n)]
        k = rng.randint(0, n + 2)
        assert top_k(a, k, largest=True) == sorted(a, reverse=True)[:min(k, n)]
        assert top_k(a, k, largest=False) == sorted(a)[:min(k, n)]


def test_known():
    assert kth_smallest([5, 2, 8, 1, 9], 1) == 2
    assert median([1, 2, 3, 4]) == 2.5
    assert top_k([5, 2, 8, 1, 9, 3], 3) == [9, 8, 5]


def test_validation():
    with pytest.raises(ValueError):
        kth_smallest([], 0)
    with pytest.raises(ValueError):
        kth_smallest([1, 2], 5)
    with pytest.raises(ValueError):
        top_k([1, 2], -1)
