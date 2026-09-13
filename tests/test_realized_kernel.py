"""Realized-kernel integrated-variance estimator."""

import math
import random

import pytest

from quantforge import realized_kernel
from quantforge.two_scale_rv import realized_variance_naive


def _bm_path(sigma, n, dt, noise, seed):
    rng = random.Random(seed)
    p = [0.0]
    x = 0.0
    for _ in range(n):
        x += sigma * math.sqrt(dt) * rng.gauss(0, 1)
        p.append(x + (noise * rng.gauss(0, 1) if noise > 0 else 0.0))
    return p


def test_no_noise_matches_realized_variance():
    p = _bm_path(0.02, 2000, 1.0 / 2000, 0.0, 42)
    rk = realized_kernel(p)
    rv = realized_variance_naive(p)
    assert abs(rk - rv) < 0.1 * rv        # close without noise


def test_no_noise_near_true_iv():
    # Single-path realized-kernel sampling error is large; average over paths.
    sigma, n, T = 0.02, 3000, 1.0
    iv = sigma * sigma * T
    est = [realized_kernel(_bm_path(sigma, n, T / n, 0.0, s)) for s in range(50)]
    mean = sum(est) / len(est)
    assert abs(mean - iv) < 0.1 * iv


def test_noise_robustness_beats_naive():
    sigma, n, T, noise = 0.02, 2000, 1.0, 0.001
    iv = sigma * sigma * T
    naive_bias, kernel_bias = [], []
    for s in range(30):
        p = _bm_path(sigma, n, T / n, noise, s)
        naive_bias.append(realized_variance_naive(p) - iv)
        kernel_bias.append(realized_kernel(p) - iv)
    mean_naive = sum(naive_bias) / len(naive_bias)
    mean_kernel = sum(kernel_bias) / len(kernel_bias)
    assert abs(mean_kernel) < abs(mean_naive)          # kernel corrects the bias
    assert abs(mean_kernel) < 0.25 * iv                # and stays close to the truth


def test_non_negative():
    for s in range(10):
        p = _bm_path(0.02, 1000, 1.0 / 1000, 0.001, s)
        assert realized_kernel(p) >= 0.0


def test_explicit_bandwidth():
    p = _bm_path(0.02, 1000, 1.0 / 1000, 0.0005, 3)
    v = realized_kernel(p, bandwidth=20)
    assert v > 0.0


def test_validation():
    with pytest.raises(ValueError):
        realized_kernel([0.0, 0.1])       # < 3 prices
