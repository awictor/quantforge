"""Tests for interval scheduling, cross-checked against brute-force subset enumeration."""

import itertools
import random

import pytest

from quantforge.scheduling import (
    weighted_interval_schedule,
    activity_selection,
    min_rooms,
)


def _overlaps(a, b):
    return a[0] < b[1] and b[0] < a[1]


def _valid_subset(ivs):
    for i in range(len(ivs)):
        for j in range(i + 1, len(ivs)):
            if _overlaps(ivs[i], ivs[j]):
                return False
    return True


def _brute_maxweight(ivs):
    best = 0.0
    for r in range(len(ivs) + 1):
        for combo in itertools.combinations(ivs, r):
            if _valid_subset(combo):
                best = max(best, sum(c[2] for c in combo))
    return best


def _brute_maxcount(ivs):
    best = 0
    for r in range(len(ivs) + 1):
        for combo in itertools.combinations(ivs, r):
            if _valid_subset(combo):
                best = max(best, len(combo))
    return best


def _brute_rooms(ivs):
    if not ivs:
        return 0
    peak = 0
    for t in {s for s, e in ivs}:
        peak = max(peak, sum(1 for s, e in ivs if s <= t < e))
    return peak


def test_fuzz_max_weight():
    rng = random.Random(151)
    for _ in range(3000):
        ivs = []
        for _ in range(rng.randint(0, 8)):
            s = rng.randint(0, 10)
            e = s + rng.randint(1, 6)
            ivs.append((s, e, rng.randint(0, 10)))
        tw, chosen = weighted_interval_schedule(ivs)
        assert tw == _brute_maxweight(ivs)
        assert _valid_subset(chosen)
        assert abs(sum(c[2] for c in chosen) - tw) < 1e-9


def test_fuzz_activity_count():
    rng = random.Random(152)
    for _ in range(3000):
        ivs = []
        for _ in range(rng.randint(0, 8)):
            s = rng.randint(0, 10)
            e = s + rng.randint(1, 6)
            ivs.append((s, e))
        sel = activity_selection(ivs)
        assert len(sel) == _brute_maxcount([(s, e, 1) for s, e in ivs])
        assert _valid_subset(sel)


def test_fuzz_min_rooms():
    rng = random.Random(153)
    for _ in range(3000):
        ivs = []
        for _ in range(rng.randint(0, 8)):
            s = rng.randint(0, 10)
            e = s + rng.randint(1, 6)
            ivs.append((s, e))
        assert min_rooms(ivs) == _brute_rooms(ivs)


def test_weighted_explicit():
    ivs = [(1, 3, 5), (2, 5, 6), (4, 6, 5), (6, 7, 4)]
    tw, chosen = weighted_interval_schedule(ivs)
    assert tw == 14
    assert chosen == [(1, 3, 5), (4, 6, 5), (6, 7, 4)]


def test_weighted_half_open_no_conflict():
    assert weighted_interval_schedule([(0, 1, 3), (1, 2, 4)])[0] == 7


def test_weighted_prefers_single_heavy():
    # one heavy interval outweighs two light overlapping ones
    tw, chosen = weighted_interval_schedule([(0, 10, 100), (0, 5, 40), (5, 10, 40)])
    assert tw == 100
    assert chosen == [(0, 10, 100)]


def test_activity_classic():
    ivs = [(1, 4), (3, 5), (0, 6), (5, 7), (3, 9), (5, 9),
           (6, 10), (8, 11), (8, 12), (2, 14), (12, 16)]
    assert activity_selection(ivs) == [(1, 4), (5, 7), (8, 11), (12, 16)]


def test_min_rooms_known():
    assert min_rooms([(0, 30), (5, 10), (15, 20)]) == 2
    assert min_rooms([(1, 2), (2, 3), (3, 4)]) == 1  # half-open back-to-back
    assert min_rooms([(0, 5), (0, 5), (0, 5)]) == 3


def test_empty_inputs():
    assert weighted_interval_schedule([]) == (0.0, [])
    assert activity_selection([]) == []
    assert min_rooms([]) == 0


def test_single_interval():
    assert weighted_interval_schedule([(0, 5, 7)]) == (7, [(0, 5, 7)])
    assert activity_selection([(0, 5)]) == [(0, 5)]
    assert min_rooms([(0, 5)]) == 1


def test_negative_weight_raises():
    with pytest.raises(ValueError):
        weighted_interval_schedule([(0, 1, -1)])


def test_zero_weight_intervals():
    # zero weights add nothing; total is 0 (any/empty selection is optimal)
    tw, chosen = weighted_interval_schedule([(0, 1, 0), (2, 3, 0)])
    assert tw == 0
