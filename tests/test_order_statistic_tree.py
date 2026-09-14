"""Tests for OrderStatisticTree, cross-checked against a sorted multiset."""

import bisect
import random

import pytest

from quantforge.order_statistic_tree import OrderStatisticTree


def test_fuzz_vs_sorted_list():
    rng = random.Random(521)
    for _ in range(3000):
        universe = sorted(set(rng.randint(-20, 20) for _ in range(rng.randint(1, 15))))
        ost = OrderStatisticTree(universe)
        ref = []
        for _ in range(30):
            op = rng.random()
            if op < 0.5 or not ref:
                v = rng.choice(universe)
                ost.add(v)
                bisect.insort(ref, v)
            elif op < 0.7:
                v = rng.choice(ref)
                ost.remove(v)
                ref.pop(bisect.bisect_left(ref, v))
            assert len(ost) == len(ref)
            x = rng.randint(-25, 25)
            assert ost.count_less(x) == bisect.bisect_left(ref, x)
            lo = rng.randint(-25, 10)
            hi = lo + rng.randint(0, 20)
            assert ost.count_range(lo, hi) == bisect.bisect_left(ref, hi) - bisect.bisect_left(ref, lo)
            cv = rng.choice(universe)
            assert ost.count(cv) == ref.count(cv)
            if ref:
                k = rng.randint(0, len(ref) - 1)
                assert ost.select(k) == ref[k]


def test_explicit():
    ost = OrderStatisticTree(range(10))
    for v in [3, 1, 4, 1, 5, 9, 2, 6]:
        ost.add(v)
    assert len(ost) == 8
    assert ost.count_less(4) == 4          # 1, 1, 2, 3
    assert ost.rank(5) == 5                # 1, 1, 2, 3, 4
    assert [ost.select(i) for i in range(8)] == [1, 1, 2, 3, 4, 5, 6, 9]
    assert ost.count(1) == 2
    assert ost.count_range(2, 6) == 4      # 2, 3, 4, 5


def test_add_remove_multiplicity():
    ost = OrderStatisticTree(range(5))
    ost.add(2, 3)
    assert ost.count(2) == 3
    assert len(ost) == 3
    ost.remove(2, 2)
    assert ost.count(2) == 1


def test_select_all_orders():
    ost = OrderStatisticTree(range(20))
    vals = [5, 5, 1, 9, 3, 3, 3, 7]
    for v in vals:
        ost.add(v)
    srt = sorted(vals)
    assert [ost.select(i) for i in range(len(srt))] == srt


def test_count_range_half_open():
    ost = OrderStatisticTree(range(10))
    for v in range(10):
        ost.add(v)
    assert ost.count_range(2, 5) == 3  # 2, 3, 4 (5 excluded)
    assert ost.count_range(0, 10) == 10
    assert ost.count_range(5, 5) == 0


def test_count_of_absent_value():
    ost = OrderStatisticTree(range(5))
    ost.add(2)
    assert ost.count(3) == 0
    assert ost.count(99) == 0  # outside universe


def test_empty_tree():
    ost = OrderStatisticTree([1, 2, 3])
    assert len(ost) == 0
    assert ost.count_less(2) == 0
    assert ost.count_range(1, 3) == 0


def test_value_outside_universe_raises():
    ost = OrderStatisticTree([1, 2, 3])
    with pytest.raises(ValueError):
        ost.add(99)


def test_remove_more_than_present_raises():
    ost = OrderStatisticTree([1, 2, 3])
    ost.add(2)
    with pytest.raises(ValueError):
        ost.remove(2, 5)


def test_select_out_of_range_raises():
    ost = OrderStatisticTree([1, 2, 3])
    ost.add(1)
    with pytest.raises(IndexError):
        ost.select(5)


def test_bad_range_raises():
    ost = OrderStatisticTree([1, 2, 3])
    with pytest.raises(ValueError):
        ost.count_range(5, 2)


def test_float_universe():
    ost = OrderStatisticTree([0.5, 1.5, 2.5, 3.5])
    ost.add(1.5)
    ost.add(3.5)
    assert ost.count_less(2.0) == 1
    assert ost.select(1) == 3.5
