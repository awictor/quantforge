"""Tests for XorBasis (GF(2) linear basis), cross-checked against subset-XOR enumeration."""

import random

import pytest

from quantforge.xor_basis import XorBasis


def _brute_reachable(vals):
    s = {0}
    for v in vals:
        s |= {x ^ v for x in s}
    return sorted(s)


def test_fuzz_vs_brute_reachable_set():
    rng = random.Random(101)
    for _ in range(4000):
        n = rng.randint(0, 8)
        vals = [rng.randint(0, 63) for _ in range(n)]
        xb = XorBasis(vals)
        reach = _brute_reachable(vals)
        reachset = set(reach)
        assert xb.count_distinct() == len(reach)
        assert xb.max_xor() == max(reach)
        assert xb.min_xor() == min(reach)
        assert (1 << xb.rank()) == len(reach)
        for t in range(64):
            assert xb.can_represent(t) == (t in reachset)
        for k in range(len(reach)):
            assert xb.kth_smallest(k) == reach[k]


def test_fuzz_with_start_offset():
    rng = random.Random(102)
    for _ in range(2000):
        vals = [rng.randint(0, 63) for _ in range(rng.randint(0, 6))]
        xb = XorBasis(vals)
        reach = _brute_reachable(vals)
        start = rng.randint(0, 63)
        assert xb.max_xor(start) == max(start ^ x for x in reach)
        assert xb.min_xor(start) == min(start ^ x for x in reach)


def test_insert_returns_independence():
    xb = XorBasis()
    assert xb.insert(5) is True
    assert xb.insert(5) is False  # already representable
    assert xb.insert(3) is True
    assert xb.insert(6) is False  # 5 ^ 3 == 6


def test_rank_and_count():
    xb = XorBasis([5, 3, 6])
    assert xb.rank() == 2
    assert xb.count_distinct() == 4
    assert len(xb) == 2


def test_max_xor_known():
    xb = XorBasis([1, 2, 4])
    assert xb.max_xor() == 7  # full 3-bit span
    assert xb.min_xor() == 0


def test_can_represent():
    xb = XorBasis([1, 2])
    assert xb.can_represent(3) is True
    assert xb.can_represent(0) is True
    assert xb.can_represent(4) is False
    assert xb.can_represent(-1) is False


def test_kth_smallest_order():
    xb = XorBasis([1, 2, 4])
    assert [xb.kth_smallest(k) for k in range(8)] == list(range(8))


def test_merge():
    a = XorBasis([1, 2])
    b = XorBasis([4])
    a.merge(b)
    assert a.max_xor() == 7
    assert a.count_distinct() == 8


def test_empty_basis():
    e = XorBasis()
    assert e.max_xor() == 0
    assert e.min_xor() == 0
    assert e.count_distinct() == 1
    assert e.rank() == 0
    assert e.kth_smallest(0) == 0
    assert e.can_represent(0) is True
    assert e.can_represent(1) is False


def test_duplicate_and_zero_inserts():
    xb = XorBasis()
    assert xb.insert(0) is False  # zero never enlarges the basis
    assert xb.rank() == 0
    xb.insert(7)
    assert xb.insert(7) is False


def test_negative_insert_raises():
    with pytest.raises(ValueError):
        XorBasis().insert(-1)


def test_kth_out_of_range_raises():
    xb = XorBasis([1])
    with pytest.raises(IndexError):
        xb.kth_smallest(2)  # only {0, 1} reachable
    with pytest.raises(IndexError):
        xb.kth_smallest(-1)


def test_merge_wrong_type_raises():
    with pytest.raises(TypeError):
        XorBasis().merge([1, 2, 3])


def test_large_values():
    vals = [1 << 40, 1 << 41, (1 << 40) | (1 << 41)]
    xb = XorBasis(vals)
    assert xb.rank() == 2  # third is the XOR of the first two
    assert xb.max_xor() == (1 << 40) | (1 << 41)
