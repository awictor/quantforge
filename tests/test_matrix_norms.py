"""Condition number, matrix rank, and matrix norms."""

import random

import pytest

from quantforge import (condition_number, matrix_rank, spectral_norm,
                        frobenius_norm)
from quantforge.svd import svd


def test_condition_of_identity():
    assert abs(condition_number([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) - 1.0) < 1e-9


def test_condition_of_scaled_diagonal():
    assert abs(condition_number([[1000, 0], [0, 1]]) - 1000.0) < 1e-6


def test_condition_singular_infinite():
    assert condition_number([[1, 2], [2, 4]]) == float("inf")


def test_rank():
    assert matrix_rank([[1, 0], [0, 1]]) == 2
    assert matrix_rank([[1, 2], [2, 4]]) == 1
    assert matrix_rank([[0, 0], [0, 0]]) == 0


def test_spectral_norm_is_sigma_max():
    rng = random.Random(1)
    A = [[rng.gauss(0, 1) for _ in range(3)] for _ in range(4)]
    _, s, _ = svd(A)
    assert abs(spectral_norm(A) - s[0]) < 1e-12


def test_frobenius_equals_root_sum_squares_of_singular_values():
    rng = random.Random(2)
    A = [[rng.gauss(0, 1) for _ in range(3)] for _ in range(4)]
    _, s, _ = svd(A)
    assert abs(frobenius_norm(A) - sum(x * x for x in s) ** 0.5) < 1e-10


def test_frobenius_direct():
    assert abs(frobenius_norm([[3.0, 4.0]]) - 5.0) < 1e-12
