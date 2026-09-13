"""Multivariate normal density and log-determinant."""

import math

import random

import pytest

from quantforge import log_determinant, mvn_logpdf, mvn_pdf
from quantforge.lu import determinant


def test_log_determinant_diagonal():
    assert abs(log_determinant([[2, 0, 0], [0, 3, 0], [0, 0, 4]]) - math.log(24)) < 1e-9


def test_log_determinant_matches_lu():
    rng = random.Random(1)
    n = 5
    A = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    S = [[sum(A[i][k] * A[j][k] for k in range(n)) + (1.0 if i == j else 0.0)
          for j in range(n)] for i in range(n)]
    assert abs(log_determinant(S) - math.log(determinant(S))) < 1e-8


def test_univariate_reduction():
    sig2 = 4.0
    lp = mvn_logpdf([0.0], [0.0], [[sig2]])
    assert abs(lp - (-0.5 * math.log(2 * math.pi * sig2))) < 1e-12


def test_two_d_integrates_to_one():
    cov = [[1.0, 0.5], [0.5, 2.0]]
    mean = [0.0, 0.0]
    N, lo, hi = 300, -8.0, 8.0
    dx = (hi - lo) / N
    total = 0.0
    for i in range(N):
        for j in range(N):
            x = lo + (i + 0.5) * dx
            y = lo + (j + 0.5) * dx
            total += mvn_pdf([x, y], mean, cov) * dx * dx
    assert abs(total - 1.0) < 1e-3


def test_peak_at_mean():
    cov = [[1.0, 0.3], [0.3, 1.0]]
    assert mvn_logpdf([0, 0], [0, 0], cov) > mvn_logpdf([1, 1], [0, 0], cov)


def test_validation():
    with pytest.raises(ValueError):
        log_determinant([[1.0, 2.0], [2.0, 1.0]])       # not positive definite
    with pytest.raises(ValueError):
        mvn_logpdf([0.0, 0.0], [0.0], [[1.0]])          # dimension mismatch
