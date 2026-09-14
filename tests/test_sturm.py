"""Tests for Sturm-sequence real-root counting/isolation, cross-checked against known roots."""

import random

import pytest

from quantforge.sturm import sturm_sequence, real_root_count, isolate_real_roots


def _poly_from_roots(roots):
    c = [1.0]
    for r in roots:
        nc = [0.0] * (len(c) + 1)
        for i, v in enumerate(c):
            nc[i] += v
            nc[i + 1] -= r * v
        c = nc
    return c


def _eval(p, x):
    v = 0.0
    for c in p:
        v = v * x + c
    return v


def test_fuzz_root_count_and_isolation():
    rng = random.Random(251)
    for _ in range(3000):
        k = rng.randint(1, 5)
        roots = sorted(set(rng.randint(-8, 8) for _ in range(k)))
        p = _poly_from_roots(roots)
        a, b = -10.0, 10.0
        assert real_root_count(p, a, b) == sum(1 for r in roots if a < r <= b)
        lo = rng.randint(-9, 0)
        hi = rng.randint(1, 9)
        assert real_root_count(p, lo, hi) == sum(1 for r in roots if lo < r <= hi)
        ivs = isolate_real_roots(p, a, b)
        expected = [r for r in roots if a < r <= b]
        assert len(ivs) == len(expected)
        for lo2, hi2 in ivs:
            inside = [r for r in roots if lo2 < r <= hi2 + 1e-9]
            assert len(inside) == 1


def test_cubic_three_roots():
    p = [1, -6, 11, -6]  # (x-1)(x-2)(x-3)
    assert real_root_count(p, 0, 4) == 3
    assert real_root_count(p, 1.5, 2.5) == 1
    assert real_root_count(p, 4, 10) == 0
    assert len(isolate_real_roots(p, 0, 4)) == 3


def test_half_open_interval_semantics():
    p = [1, -6, 11, -6]
    assert real_root_count(p, 0, 1) == 1  # (0, 1] includes the root at 1
    assert real_root_count(p, 1, 2) == 1  # (1, 2] includes 2, excludes 1


def test_no_real_roots():
    assert real_root_count([1, 0, 1], -10, 10) == 0  # x^2 + 1
    assert isolate_real_roots([1, 0, 1], -10, 10) == []


def test_repeated_root_counted_once():
    assert real_root_count([1, -4, 4], 0, 5) == 1  # (x-2)^2


def test_isolation_brackets_sign_change():
    p = [1, -6, 11, -6]
    for lo, hi in isolate_real_roots(p, 0, 4):
        # a bracketing interval has opposite signs at its ends (or a root at an end)
        assert _eval(p, lo) * _eval(p, hi) <= 1e-6


def test_linear_polynomial():
    assert real_root_count([1, -3], 0, 5) == 1  # x - 3
    assert real_root_count([1, -3], 4, 5) == 0


def test_sturm_sequence_starts_with_poly_and_derivative():
    seq = sturm_sequence([1, -6, 11, -6])
    assert seq[0] == [1.0, -6.0, 11.0, -6.0]
    assert seq[1] == [3.0, -12.0, 11.0]  # derivative


def test_negative_root_interval():
    p = _poly_from_roots([-5, -2, 3])
    assert real_root_count(p, -10, 10) == 3
    assert real_root_count(p, -10, 0) == 2
    assert real_root_count(p, 0, 10) == 1


def test_a_ge_b_raises():
    with pytest.raises(ValueError):
        real_root_count([1, -1], 5, 1)
    with pytest.raises(ValueError):
        isolate_real_roots([1, -1], 5, 1)


def test_high_degree_distinct_roots():
    roots = [-7, -3, -1, 2, 4, 6]
    p = _poly_from_roots(roots)
    assert real_root_count(p, -10, 10) == 6
    assert len(isolate_real_roots(p, -10, 10)) == 6
