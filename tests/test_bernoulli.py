"""Bernoulli numbers and Faulhaber's sum-of-powers formula."""

from fractions import Fraction

import pytest

from quantforge import bernoulli_number, faulhaber, bernoulli_sequence


def test_known_bernoulli_values():
    known = {0: Fraction(1), 1: Fraction(1, 2), 2: Fraction(1, 6),
             4: Fraction(-1, 30), 6: Fraction(1, 42), 8: Fraction(-1, 30),
             10: Fraction(5, 66)}
    for k, v in known.items():
        assert bernoulli_number(k) == v
    for k in (3, 5, 7, 9, 11):
        assert bernoulli_number(k) == 0            # odd > 1 vanish


def test_faulhaber_matches_direct_sum():
    for p in range(8):
        for m in (0, 1, 2, 5, 10, 50):
            assert faulhaber(m, p) == sum(k ** p for k in range(1, m + 1))


def test_faulhaber_closed_forms():
    for m in (1, 7, 20, 100):
        assert faulhaber(m, 1) == m * (m + 1) // 2
        assert faulhaber(m, 2) == m * (m + 1) * (2 * m + 1) // 6
        assert faulhaber(m, 3) == (m * (m + 1) // 2) ** 2
    assert faulhaber(10, 2) == 385
    assert faulhaber(10, 3) == 3025


def test_sequence():
    seq = bernoulli_sequence(6)
    assert seq[0] == 1 and seq[1] == Fraction(1, 2) and seq[2] == Fraction(1, 6)
    assert len(seq) == 7


def test_validation():
    with pytest.raises(ValueError):
        bernoulli_number(-1)
    with pytest.raises(ValueError):
        faulhaber(-1, 2)
    with pytest.raises(ValueError):
        faulhaber(5, -1)
