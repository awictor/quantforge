"""Two-scale realized variance (microstructure-noise robust)."""

import random

import pytest

from quantforge import (
    two_scale_realized_variance, realized_variance_naive, noise_variance_estimate,
)


def _efficient_walk(n, sigma, seed):
    rng = random.Random(seed)
    p = [0.0]
    for _ in range(n):
        p.append(p[-1] + rng.gauss(0, sigma))
    return p


def _add_noise(p, noise_sd, seed):
    rng = random.Random(seed)
    return [pi + rng.gauss(0, noise_sd) for pi in p]


def test_noiseless_tsrv_near_integrated_variance():
    n, sigma = 2000, 0.001
    p = _efficient_walk(n, sigma, 1)
    iv = n * sigma * sigma
    assert abs(two_scale_realized_variance(p) - iv) / iv < 0.15


def test_tsrv_beats_naive_under_noise():
    n, sigma, noise_sd = 2000, 0.001, 0.0005
    p = _efficient_walk(n, sigma, 1)
    pobs = _add_noise(p, noise_sd, 2)
    iv = n * sigma * sigma
    naive = realized_variance_naive(pobs)
    ts = two_scale_realized_variance(pobs)
    assert abs(ts - iv) < abs(naive - iv)         # TSRV less biased
    assert naive > iv                              # naive biased upward


def test_naive_bias_is_two_n_noise_variance():
    n, sigma, noise_sd = 2000, 0.001, 0.0005
    p = _efficient_walk(n, sigma, 1)
    pobs = _add_noise(p, noise_sd, 2)
    iv = n * sigma * sigma
    bias = realized_variance_naive(pobs) - iv
    assert abs(bias - 2 * n * noise_sd ** 2) / (2 * n * noise_sd ** 2) < 0.2


def test_noise_variance_estimate_in_noise_dominated_regime():
    # Estimator assumes noise dominates the finest scale (real tick regime).
    n, sigma, noise_sd = 5000, 0.0002, 0.002
    p = _efficient_walk(n, sigma, 3)
    pobs = _add_noise(p, noise_sd, 4)
    est = noise_variance_estimate(pobs)
    assert abs(est - noise_sd ** 2) / noise_sd ** 2 < 0.15


def test_validation():
    with pytest.raises(ValueError):
        two_scale_realized_variance([1.0, 2.0])       # n < 3
    with pytest.raises(ValueError):
        two_scale_realized_variance([1.0, 2.0, 3.0], K=1)
    with pytest.raises(ValueError):
        realized_variance_naive([1.0])
    with pytest.raises(ValueError):
        noise_variance_estimate([1.0])
