"""Tests for matrix-power linear recurrences, cross-checked against direct iteration."""

import random

import pytest

from quantforge.linear_recurrence import matrix_power, linear_recurrence_nth, fibonacci


def _direct_rec(coeffs, initial, n, mod=None):
    k = len(coeffs)
    if n < k:
        return initial[n] % mod if mod else initial[n]
    x = list(initial)
    for i in range(k, n + 1):
        v = sum(coeffs[j] * x[i - 1 - j] for j in range(k))
        x.append(v % mod if mod else v)
    return x[n] % mod if mod else x[n]


def _mat_mul(a, b, mod=None):
    n = len(a)
    out = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = sum(a[i][t] * b[t][j] for t in range(n))
            out[i][j] = s % mod if mod else s
    return out


def _direct_matpow(m, p, mod=None):
    n = len(m)
    r = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    for _ in range(p):
        r = _mat_mul(r, m, mod)
    return r


def test_fibonacci_known_values():
    fibs = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]
    for i, f in enumerate(fibs):
        assert fibonacci(i) == f
    assert fibonacci(50) == 12586269025
    assert fibonacci(100) == 354224848179261915075


def test_fuzz_recurrence_vs_direct():
    rng = random.Random(191)
    for _ in range(2000):
        k = rng.randint(1, 4)
        coeffs = [rng.randint(-3, 3) for _ in range(k)]
        initial = [rng.randint(-5, 5) for _ in range(k)]
        n = rng.randint(0, 40)
        assert linear_recurrence_nth(coeffs, initial, n) == _direct_rec(coeffs, initial, n)
        mod = rng.choice([None, 1000, 998244353])
        assert linear_recurrence_nth(coeffs, initial, n, mod=mod) == _direct_rec(
            coeffs, initial, n, mod
        )


def test_fuzz_matrix_power_vs_repeated():
    rng = random.Random(192)
    for _ in range(1000):
        k = rng.randint(1, 4)
        m = [[rng.randint(-3, 3) for _ in range(k)] for _ in range(k)]
        p = rng.randint(0, 12)
        mod = rng.choice([None, 97])
        assert matrix_power(m, p, mod=mod) == _direct_matpow(m, p, mod)


def test_tribonacci():
    seq = [0, 0, 1, 1, 2, 4, 7, 13, 24, 44, 81]
    for i, val in enumerate(seq):
        assert linear_recurrence_nth([1, 1, 1], [0, 0, 1], i) == val


def test_fibonacci_mod():
    assert fibonacci(1000, mod=1000000007) == _direct_rec([1, 1], [0, 1], 1000, 1000000007)


def test_n_below_k_returns_initial():
    assert linear_recurrence_nth([1, 1], [5, 7], 0) == 5
    assert linear_recurrence_nth([1, 1], [5, 7], 1) == 7
    assert linear_recurrence_nth([1, 1], [5, 7], 0, mod=3) == 5 % 3


def test_matrix_power_zero_is_identity():
    assert matrix_power([[2, 3], [4, 5]], 0) == [[1, 0], [0, 1]]


def test_matrix_power_one_is_self():
    m = [[2, 3], [4, 5]]
    assert matrix_power(m, 1) == m


def test_geometric_scalar_recurrence():
    # x_n = 2 x_{n-1}, x_0 = 1 -> powers of two
    for n in range(15):
        assert linear_recurrence_nth([2], [1], n) == 2 ** n


def test_large_n_is_fast_and_correct():
    # matrix power handles a huge index that direct iteration also can (checked mod)
    mod = 1000000007
    assert fibonacci(10**6, mod=mod) == _direct_rec([1, 1], [0, 1], 10**6, mod)


def test_non_square_matrix_raises():
    with pytest.raises(ValueError):
        matrix_power([[1, 2]], 2)


def test_negative_power_raises():
    with pytest.raises(ValueError):
        matrix_power([[1]], -1)


def test_initial_length_mismatch_raises():
    with pytest.raises(ValueError):
        linear_recurrence_nth([1, 1], [1], 3)


def test_negative_n_raises():
    with pytest.raises(ValueError):
        linear_recurrence_nth([1, 1], [0, 1], -1)
    with pytest.raises(ValueError):
        fibonacci(-1)
