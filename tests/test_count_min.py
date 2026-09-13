"""Count-Min sketch and Bloom filter."""

import random

import pytest

from quantforge import CountMinSketch, BloomFilter


def test_count_min_never_underestimates():
    cms = CountMinSketch(width=2048, depth=5)
    truth = {}
    rng = random.Random(3)
    for _ in range(50000):
        k = f"item{rng.randint(0, 999)}"
        cms.add(k)
        truth[k] = truth.get(k, 0) + 1
    assert all(cms.estimate(k) >= truth[k] for k in truth)


def test_count_min_heavy_hitter_accurate():
    cms = CountMinSketch(width=4096, depth=5)
    for _ in range(10000):
        cms.add("hot")
    rng = random.Random(1)
    for _ in range(5000):
        cms.add(f"cold{rng.randint(0, 4999)}")
    assert abs(cms.estimate("hot") - 10000) / 10000 < 0.05


def test_count_min_bulk_add():
    cms = CountMinSketch()
    cms.add("x", count=100)
    assert cms.estimate("x") == 100


def test_bloom_no_false_negatives():
    bf = BloomFilter(capacity=10000, error_rate=0.01)
    members = [f"m{i}" for i in range(10000)]
    for m in members:
        bf.add(m)
    assert all(m in bf for m in members)


def test_bloom_false_positive_rate_near_target():
    bf = BloomFilter(capacity=10000, error_rate=0.01)
    for i in range(10000):
        bf.add(f"m{i}")
    fp = sum(1 for i in range(10000) if f"nonmember{i}" in bf)
    assert fp / 10000 < 0.03            # near the 0.01 target, allow slack


def test_bloom_empty():
    assert "x" not in BloomFilter(100)


def test_validation():
    with pytest.raises(ValueError):
        CountMinSketch(width=0)
    with pytest.raises(ValueError):
        BloomFilter(capacity=0)
    with pytest.raises(ValueError):
        BloomFilter(100, error_rate=1.5)
