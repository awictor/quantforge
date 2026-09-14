"""Tests for Berlekamp-Massey, cross-checked by regenerating and extending sequences."""

import random

import pytest

from quantforge.berlekamp_massey import berlekamp_massey, berlekamp_massey_next

MOD = 998244353


def _gen(coeffs, init, n, mod):
    L = len(coeffs)
    s = [x % mod for x in init[:L]]
    while len(s) < n:
        s.append(sum(coeffs[j] * s[-1 - j] for j in range(L)) % mod)
    return s


def _satisfies(coeffs, s, mod):
    L = len(coeffs)
    return all(
        sum(coeffs[j] * s[i - 1 - j] for j in range(L)) % mod == s[i]
        for i in range(L, len(s))
    )


def test_fuzz_recovers_and_predicts():
    rng = random.Random(271)
    for _ in range(3000):
        L = rng.randint(1, 6)
        coeffs = [rng.randint(0, MOD - 1) for _ in range(L)]
        while coeffs[-1] == 0:
            coeffs[-1] = rng.randint(1, MOD - 1)
        init = [rng.randint(0, MOD - 1) for _ in range(L)]
        n = rng.randint(2 * L, 2 * L + 15)
        s = _gen(coeffs, init, n, MOD)
        rec = berlekamp_massey(s, MOD)
        assert _satisfies(rec, s, MOD)
        assert len(rec) <= L  # minimal, given >= 2L terms
        true_next = _gen(coeffs, init, n + 5, MOD)[n:]
        assert berlekamp_massey_next(s, MOD, 5) == true_next


def test_fibonacci():
    fib = [0, 1]
    for _ in range(20):
        fib.append(fib[-1] + fib[-2])
    assert berlekamp_massey(fib, MOD) == [1, 1]
    assert berlekamp_massey_next(fib, MOD, 1)[0] == (fib[-1] + fib[-2]) % MOD


def test_geometric():
    assert berlekamp_massey([1, 2, 4, 8, 16], MOD) == [2]
    assert berlekamp_massey_next([1, 2, 4, 8, 16], MOD, 2) == [32, 64]


def test_constant():
    assert berlekamp_massey([5, 5, 5, 5], MOD) == [1]


def test_all_zeros():
    assert berlekamp_massey([0, 0, 0, 0], MOD) == []
    assert berlekamp_massey_next([0, 0, 0], MOD, 3) == [0, 0, 0]


def test_tribonacci():
    trib = [0, 0, 1]
    for _ in range(15):
        trib.append(trib[-1] + trib[-2] + trib[-3])
    assert berlekamp_massey(trib, MOD) == [1, 1, 1]


def test_interleaved_zero_sequence():
    # a sequence with zeros in odd positions stressed the minimality bug that was fixed
    rng = random.Random(99)
    coeffs = [0, rng.randint(1, MOD - 1)]  # s_i = c * s_{i-2}
    init = [rng.randint(1, MOD - 1), 0]
    s = _gen(coeffs, init, 20, MOD)
    rec = berlekamp_massey(s, MOD)
    assert _satisfies(rec, s, MOD)
    assert len(rec) <= 2


def test_predict_matches_manual_roll():
    s = [1, 1, 2, 3, 5, 8, 13]
    coeffs = berlekamp_massey(s, MOD)
    nxt = berlekamp_massey_next(s, MOD, 3)
    assert nxt == [21, 34, 55]


def test_zero_count_prediction():
    assert berlekamp_massey_next([1, 2, 4, 8], MOD, 0) == []


def test_negative_count_raises():
    with pytest.raises(ValueError):
        berlekamp_massey_next([1, 2, 4], MOD, -1)


def test_small_prime_field():
    mod = 7
    s = _gen([3, 5], [1, 4], 12, mod)
    rec = berlekamp_massey(s, mod)
    assert _satisfies(rec, s, mod)
    assert len(rec) <= 2
