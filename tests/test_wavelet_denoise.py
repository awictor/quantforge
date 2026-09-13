"""Wavelet denoising by Haar coefficient shrinkage."""

import math
import random

import pytest

from quantforge import (
    soft_threshold,
    hard_threshold,
    mad_sigma,
    universal_threshold,
    wavelet_denoise,
    haar_dwt,
    haar_idwt,
)


def test_soft_threshold():
    assert soft_threshold(3.0, 1.0) == 2.0
    assert soft_threshold(-3.0, 1.0) == -2.0
    assert soft_threshold(0.5, 1.0) == 0.0


def test_hard_threshold():
    assert hard_threshold(3.0, 1.0) == 3.0
    assert hard_threshold(0.5, 1.0) == 0.0
    assert hard_threshold(-2.0, 1.0) == -2.0


def test_universal_threshold_formula():
    assert abs(universal_threshold(1024, 1.0) - math.sqrt(2 * math.log(1024))) < 1e-12
    assert universal_threshold(1, 5.0) == 0.0


def test_mad_sigma_recovers_noise_scale():
    # Haar transform is orthonormal, so the detail of iid noise keeps its scale.
    rng = random.Random(7)
    noise = [rng.gauss(0, 2.0) for _ in range(4096)]
    _, d = haar_dwt(noise)
    assert abs(mad_sigma(d[0]) - 2.0) < 0.15


def test_denoise_reduces_mse():
    rng = random.Random(11)
    clean = [math.sin(2 * math.pi * i / 128) for i in range(1024)]
    noisy = [clean[i] + rng.gauss(0, 0.4) for i in range(1024)]
    den = wavelet_denoise(noisy)
    mse_noisy = sum((noisy[i] - clean[i]) ** 2 for i in range(1024)) / 1024
    mse_den = sum((den[i] - clean[i]) ** 2 for i in range(1024)) / 1024
    assert mse_den < mse_noisy


def test_zero_threshold_is_identity():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(64)]
    den = wavelet_denoise(x, threshold=0.0)
    assert max(abs(den[i] - x[i]) for i in range(64)) < 1e-10


def test_huge_threshold_keeps_only_approximation():
    x = [math.sin(2 * math.pi * i / 32) for i in range(256)]
    big = wavelet_denoise(x, threshold=1e9)
    approx, details = haar_dwt(x)
    approx_only = haar_idwt(approx, [[0.0] * len(l) for l in details])
    assert max(abs(big[i] - approx_only[i]) for i in range(256)) < 1e-10


def test_hard_mode_runs():
    rng = random.Random(5)
    x = [rng.gauss(0, 1) for _ in range(64)]
    den = wavelet_denoise(x, mode="hard")
    assert len(den) == 64


def test_validation():
    with pytest.raises(ValueError):
        wavelet_denoise([1.0, 2.0, 3.0, 4.0], mode="bad")
