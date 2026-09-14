"""Tests for symmetric functions and Newton's identities, vs itertools/expansion."""

import itertools
import random
from fractions import Fraction

import pytest

from quantforge.symmetric import (
    elementary_symmetric,
    power_sums,
    power_to_elementary,
    elementary_to_power,
    poly_from_roots,
)


def _prod(xs):
    r = 1
    for x in xs:
        r *= x
    return r


def _brute_e(vals):
    n = len(vals)
    return [sum(_prod(c) for c in itertools.combinations(vals, k)) for k in range(n + 1)]


def _expand_from_roots(roots):
    c = [1]
    for r in roots:
        nc = [0] * (len(c) + 1)
        for i, v in enumerate(c):
            nc[i] += v
            nc[i + 1] -= r * v
        c = nc
    return c


def test_fuzz_all_conversions():
    rng = random.Random(261)
    for _ in range(4000):
        n = rng.randint(0, 7)
        vals = [rng.randint(-5, 5) for _ in range(n)]
        e = elementary_symmetric(vals)
        assert e == _brute_e(vals)
        assert poly_from_roots(vals) == _expand_from_roots(vals)
        ps = power_sums(vals, n)
        for k in range(n + 1):
            assert ps[k] == sum(x ** k for x in vals)
        if n > 0:
            assert power_to_elementary(ps[1:]) == [Fraction(x) for x in e]
            assert elementary_to_power(e) == ps[1:]


def test_elementary_symmetric_explicit():
    assert elementary_symmetric([1, 2, 3]) == [1, 6, 11, 6]
    assert elementary_symmetric([2, 2]) == [1, 4, 4]


def test_poly_from_roots_explicit():
    assert poly_from_roots([1, 2, 3]) == [1, -6, 11, -6]  # t^3 - 6t^2 + 11t - 6
    assert poly_from_roots([5]) == [1, -5]


def test_power_sums_explicit():
    assert power_sums([1, 2, 3], 3) == [3, 6, 14, 36]  # p0=3, p1=6, p2=14, p3=36


def test_newton_round_trip():
    vals = [2, 3, 5]
    e = elementary_symmetric(vals)
    p = elementary_to_power(e)
    assert power_to_elementary(p) == [Fraction(x) for x in e]


def test_power_to_elementary_exact_fractions():
    # power sums of a single value 3: p1=3, p2=9, p3=27 -> e = [1, 3, 0, 0]
    e = power_to_elementary([3, 9, 27])
    assert e == [Fraction(1), Fraction(3), Fraction(0), Fraction(0)]


def test_elementary_to_power_with_kmax():
    e = elementary_symmetric([1, 2, 3])
    p = elementary_to_power(e, kmax=5)
    assert p[:3] == [6, 14, 36]
    assert len(p) == 5


def test_empty_input():
    assert elementary_symmetric([]) == [1]
    assert power_sums([], 3) == [0, 0, 0, 0]
    assert poly_from_roots([]) == [1]


def test_fractional_roots():
    # roots 1/2 and 1/3: poly = t^2 - (5/6) t + 1/6
    e = elementary_symmetric([Fraction(1, 2), Fraction(1, 3)])
    assert e == [Fraction(1), Fraction(5, 6), Fraction(1, 6)]


def test_negative_kmax_raises():
    with pytest.raises(ValueError):
        power_sums([1, 2], -1)


def test_repeated_roots():
    # (t-2)^3 = t^3 - 6t^2 + 12t - 8
    assert poly_from_roots([2, 2, 2]) == [1, -6, 12, -8]
