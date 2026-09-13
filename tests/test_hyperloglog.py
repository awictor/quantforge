"""HyperLogLog distinct-count estimation."""

import math

import pytest

from quantforge import HyperLogLog


def test_accuracy_across_scales():
    for true_n in (100, 1000, 10000, 100000):
        hll = HyperLogLog(14)
        for i in range(true_n):
            hll.add(i)
        err = abs(hll.count() - true_n) / true_n
        assert err < 0.05          # comfortably within a few standard errors


def test_duplicates_do_not_inflate():
    hll = HyperLogLog(12)
    for _ in range(5):
        for i in range(1000):
            hll.add(i)
    assert abs(hll.count() - 1000) / 1000 < 0.05


def test_merge_is_union():
    a = HyperLogLog(12)
    b = HyperLogLog(12)
    for i in range(5000):
        a.add(i)
    for i in range(2500, 7500):
        b.add(i)
    a.merge(b)                     # true union 0..7500
    assert abs(a.count() - 7500) / 7500 < 0.05


def test_small_range_linear_counting():
    hll = HyperLogLog(14)
    for i in range(50):
        hll.add(i)
    assert abs(hll.count() - 50) < 5


def test_deterministic():
    h1 = HyperLogLog(10)
    h2 = HyperLogLog(10)
    for i in range(2000):
        h1.add(f"item{i}")
        h2.add(f"item{i}")
    assert h1.count() == h2.count()


def test_strings():
    hll = HyperLogLog(12)
    for i in range(3000):
        hll.add(f"user_{i}@example.com")
    assert abs(hll.count() - 3000) / 3000 < 0.05


def test_validation():
    with pytest.raises(ValueError):
        HyperLogLog(3)
    with pytest.raises(ValueError):
        HyperLogLog(12).merge(HyperLogLog(10))
