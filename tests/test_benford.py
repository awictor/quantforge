"""Benford's law first-digit analysis."""

import random

import pytest

from quantforge import (benford_expected, first_digit, first_digit_distribution,
                        benford_chi_square, benford_mad)


def test_expected_sums_to_one():
    exp = benford_expected()
    assert abs(sum(exp) - 1.0) < 1e-12
    assert abs(exp[0] - 0.30103) < 1e-4       # P(1)
    assert abs(exp[8] - 0.04576) < 1e-4       # P(9)


def test_first_digit():
    assert first_digit(0.0037) == 3
    assert first_digit(12345) == 1
    assert first_digit(9.9) == 9
    assert first_digit(-800) == 8


def test_distribution_sums():
    counts, prop = first_digit_distribution([1, 2, 2, 30, 400])
    assert sum(counts) == 5
    assert abs(sum(prop) - 1.0) < 1e-12


def test_fibonacci_conforms():
    fib = [1, 1]
    for _ in range(500):
        fib.append(fib[-1] + fib[-2])
    _, p = benford_chi_square(fib)
    assert p > 0.05
    assert benford_mad(fib) < 0.006          # Nigrini "close conformance"


def test_powers_of_two_conform():
    pw = [2 ** k for k in range(1, 600)]
    _, p = benford_chi_square(pw)
    assert p > 0.05


def test_uniform_rejected():
    rng = random.Random(1)
    unif = [rng.uniform(1, 10) for _ in range(5000)]
    chi2, p = benford_chi_square(unif)
    assert p < 0.001
    assert benford_mad(unif) > 0.015          # nonconformity


def test_validation():
    with pytest.raises(ValueError):
        first_digit(0.0)
    with pytest.raises(ValueError):
        first_digit_distribution([0, 0, 0])
