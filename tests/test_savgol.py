"""Savitzky-Golay smoothing and differentiation filter."""

import random
import statistics

import pytest

from quantforge import savgol_coeffs, savgol_filter


def test_smooth_coeffs_sum_to_one():
    assert abs(sum(savgol_coeffs(5, 2, 0)) - 1.0) < 1e-12
    assert abs(sum(savgol_coeffs(7, 3, 0)) - 1.0) < 1e-12


def test_derivative_coeffs_sum_to_zero():
    assert abs(sum(savgol_coeffs(5, 2, 1))) < 1e-12


def test_reproduces_polynomial():
    poly = [2 - 3 * t + 0.5 * t * t for t in range(20)]
    sm = savgol_filter(poly, 5, 2, 0)
    assert max(abs(sm[i] - poly[i]) for i in range(20)) < 1e-9


def test_first_derivative_of_quadratic():
    poly = [2 - 3 * t + 0.5 * t * t for t in range(20)]
    d = savgol_filter(poly, 5, 2, 1)
    true_d = [-3 + t for t in range(20)]
    assert max(abs(d[i] - true_d[i]) for i in range(20)) < 1e-9


def test_reduces_noise_variance():
    rng = random.Random(1)
    noisy = [5 + rng.gauss(0, 1) for _ in range(200)]
    sm = savgol_filter(noisy, 11, 3, 0)
    assert statistics.pvariance(sm) < 0.5 * statistics.pvariance(noisy)


def test_output_length():
    data = list(range(30))
    assert len(savgol_filter(data, 7, 2, 0)) == 30


def test_validation():
    with pytest.raises(ValueError):
        savgol_coeffs(4, 2)              # even window
    with pytest.raises(ValueError):
        savgol_coeffs(5, 5)              # degree >= window
    with pytest.raises(ValueError):
        savgol_filter([1.0, 2.0], 5, 2)  # data shorter than window
