"""Tests for Stern-Brocot tree and Farey sequences, cross-checked against brute references."""

import math
import random
from fractions import Fraction
from math import gcd

import pytest

from quantforge.stern_brocot import (
    mediant,
    stern_brocot_path,
    stern_brocot_from_path,
    best_rational_bounded,
    farey_sequence,
)


def _brute_farey(n):
    s = set()
    for d in range(1, n + 1):
        for num in range(0, d + 1):
            if gcd(num, d) == 1:
                s.add(Fraction(num, d))
    return sorted(s)


def _brute_best(x, D):
    best = None
    bd = float("inf")
    for den in range(1, D + 1):
        num = round(x * den)
        for nn in (num - 1, num, num + 1):
            f = Fraction(nn, den)
            d = abs(float(f) - x)
            if d < bd - 1e-15:
                bd = d
                best = f
    return best


def test_fuzz_path_roundtrip():
    rng = random.Random(291)
    for _ in range(5000):
        t = Fraction(rng.randint(1, 50), rng.randint(1, 50))
        assert stern_brocot_from_path(stern_brocot_path(t)) == t


def test_fuzz_best_rational_bounded():
    rng = random.Random(292)
    for _ in range(2000):
        x = rng.uniform(0, 10)
        D = rng.randint(1, 40)
        got = best_rational_bounded(x, D)
        assert got.denominator <= D
        assert abs(float(got) - x) <= abs(float(_brute_best(x, D)) - x) + 1e-12


def test_fuzz_farey_vs_brute():
    for n in range(1, 15):
        assert farey_sequence(n) == _brute_farey(n)


def test_mediant():
    assert mediant((1, 2), (1, 3)) == Fraction(2, 5)
    assert mediant(Fraction(0, 1), Fraction(1, 1)) == Fraction(1, 2)
    assert mediant(Fraction(1, 3), Fraction(1, 2)) == Fraction(2, 5)


def test_path_explicit():
    assert stern_brocot_path(Fraction(1, 1)) == ""
    assert stern_brocot_path(Fraction(1, 2)) == "L"
    assert stern_brocot_path(Fraction(2, 1)) == "R"
    assert stern_brocot_path(Fraction(2, 3)) == "LR"


def test_from_path_explicit():
    assert stern_brocot_from_path("") == Fraction(1, 1)
    assert stern_brocot_from_path("L") == Fraction(1, 2)
    assert stern_brocot_from_path("RR") == Fraction(3, 1)


def test_best_rational_pi():
    assert best_rational_bounded(math.pi, 7) == Fraction(22, 7)
    assert best_rational_bounded(math.pi, 113) == Fraction(355, 113)


def test_best_rational_exact():
    assert best_rational_bounded(0.5, 10) == Fraction(1, 2)
    assert best_rational_bounded(3.0, 5) == Fraction(3, 1)


def test_farey_explicit():
    assert farey_sequence(1) == [Fraction(0, 1), Fraction(1, 1)]
    assert farey_sequence(2) == [Fraction(0, 1), Fraction(1, 2), Fraction(1, 1)]
    assert farey_sequence(3) == [
        Fraction(0, 1), Fraction(1, 3), Fraction(1, 2), Fraction(2, 3), Fraction(1, 1)
    ]


def test_farey_length():
    # |F_n| = 1 + sum_{k=1}^n phi(k); for n=5 that is 11
    assert len(farey_sequence(5)) == 11


def test_negative_target_raises():
    with pytest.raises(ValueError):
        stern_brocot_path(Fraction(-1, 2))


def test_farey_below_one_raises():
    with pytest.raises(ValueError):
        farey_sequence(0)


def test_best_rational_bad_denominator_raises():
    with pytest.raises(ValueError):
        best_rational_bounded(1.5, 0)


def test_from_path_invalid_char_raises():
    with pytest.raises(ValueError):
        stern_brocot_from_path("LXR")
