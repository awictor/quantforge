"""Tests for inversion counting and Kendall-tau distance, cross-checked against brute force."""

import itertools
import random

import pytest

from quantforge.inversions import count_inversions, kendall_tau_distance, is_sorted


def _brute_inv(a):
    return sum(1 for i in range(len(a)) for j in range(i + 1, len(a)) if a[i] > a[j])


def _brute_kendall(a, b):
    pos_b = {x: i for i, x in enumerate(b)}
    return sum(1 for x, y in itertools.combinations(a, 2) if pos_b[x] > pos_b[y])


def test_fuzz_inversions_vs_brute():
    rng = random.Random(241)
    for _ in range(5000):
        n = rng.randint(0, 30)
        a = [rng.randint(0, 15) for _ in range(n)]  # duplicates allowed
        assert count_inversions(a) == _brute_inv(a)


def test_fuzz_kendall_vs_brute():
    rng = random.Random(242)
    for _ in range(3000):
        n = rng.randint(0, 10)
        a = rng.sample(range(n), n)
        b = rng.sample(range(n), n)
        assert kendall_tau_distance(a, b) == _brute_kendall(a, b)


def test_inversions_sorted_and_reversed():
    assert count_inversions([1, 2, 3, 4, 5]) == 0
    assert count_inversions([5, 4, 3, 2, 1]) == 10  # n(n-1)/2
    assert count_inversions([]) == 0
    assert count_inversions([1]) == 0


def test_inversions_with_duplicates():
    assert count_inversions([2, 1, 3, 1, 2]) == _brute_inv([2, 1, 3, 1, 2])
    assert count_inversions([1, 1, 1]) == 0  # equal elements are not inversions


def test_kendall_identical_and_reversed():
    assert kendall_tau_distance(["a", "b", "c"], ["a", "b", "c"]) == 0
    assert kendall_tau_distance(["a", "b", "c"], ["c", "b", "a"]) == 3  # all 3 pairs flip


def test_kendall_single_swap():
    assert kendall_tau_distance([1, 2, 3, 4], [2, 1, 3, 4]) == 1


def test_kendall_max_distance():
    n = 6
    a = list(range(n))
    b = list(reversed(a))
    assert kendall_tau_distance(a, b) == n * (n - 1) // 2


def test_is_sorted():
    assert is_sorted([1, 2, 2, 3]) is True
    assert is_sorted([1, 3, 2]) is False
    assert is_sorted([1, 2, 3], strict=True) is True
    assert is_sorted([1, 2, 2], strict=True) is False
    assert is_sorted([]) is True
    assert is_sorted([5]) is True


def test_kendall_length_mismatch_raises():
    with pytest.raises(ValueError):
        kendall_tau_distance([1, 2], [1, 2, 3])


def test_kendall_different_items_raises():
    with pytest.raises(ValueError):
        kendall_tau_distance([1, 2], [1, 3])


def test_kendall_duplicates_raise():
    with pytest.raises(ValueError):
        kendall_tau_distance([1, 1], [1, 1])


def test_large_reversed_is_fast():
    # merge sort keeps this O(n log n); a brute count would be 5e9 ops
    n = 100000
    assert count_inversions(list(range(n, 0, -1))) == n * (n - 1) // 2
