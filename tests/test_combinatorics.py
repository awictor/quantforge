"""Exact integer combinatorics."""

import math

import pytest

from quantforge import (
    binomial,
    multinomial,
    stirling_second,
    bell,
    catalan,
    partition_count,
    derangements,
)


def test_binomial_sum_identity():
    for n in range(25):
        assert sum(binomial(n, k) for k in range(n + 1)) == 2 ** n
    assert binomial(5, 7) == 0
    assert binomial(100, 50) == math.comb(100, 50)


def test_multinomial():
    assert multinomial([1, 4, 4, 2]) == 34650      # "mississippi"
    assert multinomial([2, 2]) == 6
    assert multinomial([5]) == 1


def test_stirling_second_and_bell():
    assert stirling_second(4, 2) == 7
    assert stirling_second(0, 0) == 1
    assert stirling_second(3, 5) == 0
    for n in range(15):
        assert bell(n) == sum(stirling_second(n, k) for k in range(n + 1))
    assert bell(5) == 52


def test_catalan_sequence():
    assert [catalan(n) for n in range(7)] == [1, 1, 2, 5, 14, 42, 132]


def test_partition_sequence():
    assert [partition_count(n) for n in range(9)] == [1, 1, 2, 3, 5, 7, 11, 15, 22]


def test_derangements_match_closed_form():
    for n in range(1, 15):
        assert derangements(n) == round(math.factorial(n) / math.e)
    assert derangements(0) == 1


def test_validation():
    with pytest.raises(ValueError):
        binomial(-1, 0)
    with pytest.raises(ValueError):
        multinomial([1, -2])
    for fn in (bell, catalan, partition_count, derangements):
        with pytest.raises(ValueError):
            fn(-1)
    with pytest.raises(ValueError):
        stirling_second(-1, 0)
