"""Tests for the Walsh-Hadamard transform and bitwise convolutions, vs brute force."""

import random

import pytest

from quantforge.walsh_hadamard import (
    fwht,
    ifwht,
    xor_convolve,
    and_convolve,
    or_convolve,
)


def _brute(a, b, op):
    n = len(a)
    c = [0] * n
    for i in range(n):
        for j in range(n):
            c[op(i, j)] += a[i] * b[j]
    return c


def test_fuzz_xor_convolution():
    rng = random.Random(161)
    for _ in range(4000):
        n = 1 << rng.randint(0, 5)
        a = [rng.randint(-5, 5) for _ in range(n)]
        b = [rng.randint(-5, 5) for _ in range(n)]
        assert xor_convolve(a, b) == _brute(a, b, lambda i, j: i ^ j)


def test_fuzz_and_convolution():
    rng = random.Random(162)
    for _ in range(4000):
        n = 1 << rng.randint(0, 5)
        a = [rng.randint(-5, 5) for _ in range(n)]
        b = [rng.randint(-5, 5) for _ in range(n)]
        assert and_convolve(a, b) == _brute(a, b, lambda i, j: i & j)


def test_fuzz_or_convolution():
    rng = random.Random(163)
    for _ in range(4000):
        n = 1 << rng.randint(0, 5)
        a = [rng.randint(-5, 5) for _ in range(n)]
        b = [rng.randint(-5, 5) for _ in range(n)]
        assert or_convolve(a, b) == _brute(a, b, lambda i, j: i | j)


def test_fwht_roundtrip():
    rng = random.Random(164)
    for _ in range(1000):
        n = 1 << rng.randint(0, 6)
        a = [rng.randint(-10, 10) for _ in range(n)]
        rt = ifwht(fwht(a))
        assert all(abs(x - y) < 1e-9 for x, y in zip(rt, a))


def test_fwht_of_impulse_is_flat():
    assert fwht([1, 0, 0, 0]) == [1, 1, 1, 1]
    assert fwht([1, 0]) == [1, 1]


def test_xor_identity():
    # convolving with the [1, 0, ...] impulse leaves the signal unchanged
    assert xor_convolve([1, 2, 3, 4], [1, 0, 0, 0]) == [1, 2, 3, 4]


def test_explicit_values():
    a = [1, 2, 3, 4]
    b = [5, 6, 7, 8]
    assert xor_convolve(a, b) == _brute(a, b, lambda i, j: i ^ j)
    assert and_convolve(a, b) == _brute(a, b, lambda i, j: i & j)
    assert or_convolve(a, b) == _brute(a, b, lambda i, j: i | j)


def test_integer_output_type():
    # integer inputs give exact integer XOR convolution (no float division)
    out = xor_convolve([1, 2, 3, 4], [5, 6, 7, 8])
    assert all(isinstance(v, int) for v in out)


def test_single_element():
    assert xor_convolve([3], [4]) == [12]
    assert and_convolve([3], [4]) == [12]
    assert or_convolve([3], [4]) == [12]
    assert fwht([7]) == [7]


def test_non_power_of_two_raises():
    with pytest.raises(ValueError):
        fwht([1, 2, 3])
    with pytest.raises(ValueError):
        xor_convolve([1, 2, 3], [1, 2, 3])


def test_empty_raises():
    with pytest.raises(ValueError):
        fwht([])


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        xor_convolve([1, 2], [1, 2, 3, 4])
    with pytest.raises(ValueError):
        and_convolve([1, 2], [1, 2, 3, 4])
    with pytest.raises(ValueError):
        or_convolve([1, 2], [1, 2, 3, 4])
