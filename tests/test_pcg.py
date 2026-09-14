"""Statistically strong PRNGs: PCG32 and xorshift128+."""

import statistics
from collections import Counter

import pytest

from quantforge import PCG32, Xorshift128Plus


def test_pcg32_reference_vector():
    p = PCG32(seed=42, seq=54)
    out = [p.next_uint32() for _ in range(6)]
    assert out == [0xa15c02b7, 0x7b47f409, 0xba1d3330,
                   0x83d2f293, 0xbfa4784b, 0xcbed606e]


def test_pcg32_uniformity():
    xs = [PCG32(seed=1).random() for _ in range(1)]  # warm the instance
    p = PCG32(seed=1)
    xs = [p.random() for _ in range(200000)]
    assert abs(statistics.mean(xs) - 0.5) < 0.005
    assert abs(statistics.pvariance(xs) - 1 / 12) < 0.002


def test_pcg32_randint_unbiased():
    p = PCG32(seed=7)
    c = Counter(p.randint(0, 5) for _ in range(600000))
    assert all(abs(c[k] / 600000 - 1 / 6) < 0.005 for k in range(6))


def test_reproducible():
    assert [PCG32(seed=9).next_uint32() for _ in range(10)] == \
        [PCG32(seed=9).next_uint32() for _ in range(10)]
    assert [Xorshift128Plus(seed=5).next_uint64() for _ in range(10)] == \
        [Xorshift128Plus(seed=5).next_uint64() for _ in range(10)]


def test_xorshift_uniformity():
    x = Xorshift128Plus(seed=3)
    xs = [x.random() for _ in range(200000)]
    assert abs(statistics.mean(xs) - 0.5) < 0.005


def test_validation():
    with pytest.raises(ValueError):
        PCG32(1).randint(5, 3)
