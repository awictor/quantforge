"""Operations on sets of 1-D intervals."""

import random

import pytest

from quantforge import (
    merge_intervals,
    total_covered_length,
    intervals_intersection,
    intervals_union,
    max_overlap,
)


def _points(ivs):
    s = set()
    for a, b in ivs:
        s.update(range(a, b + 1))
    return s


def test_merge_known():
    assert merge_intervals([(1, 3), (2, 6), (8, 10), (15, 18)]) == [(1, 6), (8, 10), (15, 18)]
    assert merge_intervals([(1, 2), (2, 3)]) == [(1, 3)]
    assert merge_intervals([]) == []


def test_intersection_union_vs_pointset():
    rng = random.Random(1)
    for _ in range(500):
        a = [(x, x + rng.randint(0, 10)) for x in
             [rng.randint(0, 40) for _ in range(rng.randint(0, 5))]]
        b = [(x, x + rng.randint(0, 10)) for x in
             [rng.randint(0, 40) for _ in range(rng.randint(0, 5))]]
        assert _points(intervals_intersection(a, b)) == _points(a) & _points(b)
        assert _points(intervals_union(a, b)) == _points(a) | _points(b)


def test_total_covered_length():
    assert total_covered_length([(1, 4), (2, 6), (8, 10)]) == 7


def test_max_overlap():
    assert max_overlap([(1, 5), (2, 6), (4, 8), (10, 12)]) == 3

    def brute(ivs):
        best = 0
        for x in sorted({c for a, b in ivs for c in (a, b)}):
            best = max(best, sum(1 for a, b in ivs if a <= x <= b))
        return best

    rng = random.Random(2)
    for _ in range(500):
        ivs = [(x, x + rng.randint(0, 8)) for x in
               [rng.randint(0, 30) for _ in range(rng.randint(1, 8))]]
        assert max_overlap(ivs) == brute(ivs)


def test_validation():
    with pytest.raises(ValueError):
        merge_intervals([(5, 3)])
