"""Continued-fraction expansion and best rational approximation."""

import math
import random
from fractions import Fraction

import pytest

from quantforge import cf_expansion, convergents, best_rational


def test_known_expansions():
    assert cf_expansion(415 / 93)[:4] == [4, 2, 6, 7]
    assert cf_expansion(math.pi, max_terms=5) == [3, 7, 15, 1, 292]
    assert cf_expansion(math.sqrt(2), max_terms=6) == [1, 2, 2, 2, 2, 2]


def test_pi_convergents():
    conv = convergents(cf_expansion(math.pi, max_terms=6))
    assert conv[:4] == [(3, 1), (22, 7), (333, 106), (355, 113)]


def test_best_rational_matches_fraction_limit_denominator():
    rng = random.Random(1)
    for _ in range(3000):
        x = rng.uniform(-5.0, 5.0)
        D = rng.choice([10, 100, 1000, 10000, 100000])
        p, q = best_rational(x, D)
        ref = Fraction(x).limit_denominator(D)
        assert q <= D and q >= 1
        # same fraction, or an equally-close tie
        if (p, q) != (ref.numerator, ref.denominator):
            e_mine = abs(p / q - x)
            e_ref = abs(ref.numerator / ref.denominator - x)
            assert e_mine <= e_ref + 1e-18


def test_pi_best_rationals():
    assert best_rational(math.pi, 1000) == (355, 113)
    assert best_rational(math.pi, 50) == (22, 7)


def test_exact_rationals_roundtrip():
    assert best_rational(0.75, 1000) == (3, 4)
    assert best_rational(0.375, 1000) == (3, 8)
    assert best_rational(2.0, 1000) == (2, 1)


def test_convergents_are_best_in_class():
    # each convergent beats every fraction with a strictly smaller denominator
    x = math.e
    conv = convergents(cf_expansion(x, max_terms=8))
    for p, q in conv[2:6]:
        err = abs(p / q - x)
        for qq in range(1, q):
            pp = round(x * qq)
            assert abs(pp / qq - x) >= err - 1e-15


def test_validation():
    with pytest.raises(ValueError):
        convergents([])
    with pytest.raises(ValueError):
        best_rational(1.5, 0)
