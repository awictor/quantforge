"""Tests for the number-theoretic transform, cross-checked against direct integer convolution."""

import random

import pytest

from quantforge.ntt import ntt, intt, convolve_mod, NTT_PRIME


def _direct_conv(a, b, mod):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = (r[i + j] + x * y) % mod
    return r


def test_fuzz_convolve_vs_direct():
    rng = random.Random(141)
    mod = NTT_PRIME
    for _ in range(3000):
        a = [rng.randint(0, mod - 1) for _ in range(rng.randint(1, 30))]
        b = [rng.randint(0, mod - 1) for _ in range(rng.randint(1, 30))]
        assert convolve_mod(a, b) == _direct_conv(a, b, mod)


def test_ntt_roundtrip():
    rng = random.Random(142)
    mod = NTT_PRIME
    for p in (1, 2, 4, 8, 16, 32, 64):
        x = [rng.randint(0, mod - 1) for _ in range(p)]
        assert intt(ntt(x)) == x


def test_small_polynomial_product():
    # (1 + 2x + 3x^2)(4 + 5x + 6x^2) = 4 + 13x + 28x^2 + 27x^3 + 18x^4
    assert convolve_mod([1, 2, 3], [4, 5, 6]) == [4, 13, 28, 27, 18]


def test_polynomial_squaring():
    assert convolve_mod([1, 1], [1, 1]) == [1, 2, 1]
    assert convolve_mod([1, 1, 1, 1], [1, 1, 1, 1]) == [1, 2, 3, 4, 3, 2, 1]


def test_big_values_wrap_mod_prime():
    mod = NTT_PRIME
    assert convolve_mod([mod - 1], [mod - 1]) == [(mod - 1) * (mod - 1) % mod]


def test_convolution_length():
    a = [1, 2, 3, 4]
    b = [5, 6]
    assert len(convolve_mod(a, b)) == len(a) + len(b) - 1


def test_single_element_convolution():
    assert convolve_mod([5], [7]) == [35]


def test_empty_convolution():
    assert convolve_mod([], [1, 2]) == []
    assert convolve_mod([1, 2], []) == []


def test_non_power_of_two_transform_raises():
    with pytest.raises(ValueError):
        ntt([1, 2, 3])
    with pytest.raises(ValueError):
        intt([1, 2, 3])


def test_transform_reduces_input_mod():
    mod = NTT_PRIME
    # inputs above the modulus are reduced first, so the round-trip lands in [0, mod)
    x = [mod + 1, 2 * mod + 3, 4, 5]
    rt = intt(ntt(x))
    assert rt == [v % mod for v in x]


def test_identity_convolution():
    # convolving with [1] leaves a sequence unchanged
    a = [3, 1, 4, 1, 5]
    assert convolve_mod(a, [1]) == a


def test_linearity_via_shift():
    # convolving with [0, 1] shifts coefficients up by one degree
    a = [2, 3, 5]
    assert convolve_mod(a, [0, 1]) == [0, 2, 3, 5]
