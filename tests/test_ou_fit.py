"""Ornstein-Uhlenbeck parameter estimation."""

import math
import random

import pytest

from quantforge import fit_ornstein_uhlenbeck


def _simulate_ou(kappa, theta, sigma, dt, n, seed):
    b = math.exp(-kappa * dt)
    a = theta * (1 - b)
    sd = math.sqrt(sigma * sigma * (1 - b * b) / (2 * kappa))
    rng = random.Random(seed)
    x = [theta]
    for _ in range(n):
        x.append(a + b * x[-1] + rng.gauss(0, sd))
    return x


def test_recovers_known_parameters():
    kappa, theta, sigma, dt = 2.0, 5.0, 0.5, 1.0 / 52
    x = _simulate_ou(kappa, theta, sigma, dt, 20000, 1)
    p = fit_ornstein_uhlenbeck(x, dt=dt)
    assert abs(p["kappa"] - kappa) / kappa < 0.15
    assert abs(p["theta"] - theta) < 0.2
    assert abs(p["sigma"] - sigma) / sigma < 0.1


def test_half_life_relation():
    x = _simulate_ou(2.0, 5.0, 0.5, 1.0 / 52, 5000, 2)
    p = fit_ornstein_uhlenbeck(x, dt=1.0 / 52)
    assert abs(p["half_life"] - math.log(2.0) / p["kappa"]) < 1e-12


def test_faster_reversion_shorter_half_life():
    slow = fit_ornstein_uhlenbeck(_simulate_ou(1.0, 0.0, 0.3, 1.0, 20000, 3))
    fast = fit_ornstein_uhlenbeck(_simulate_ou(4.0, 0.0, 0.3, 1.0, 20000, 4))
    assert fast["half_life"] < slow["half_life"]


def test_random_walk_gives_tiny_kappa():
    rng = random.Random(1)
    rw = [0.0]
    for _ in range(2000):
        rw.append(rw[-1] + rng.gauss(0, 1))
    p = fit_ornstein_uhlenbeck(rw)
    assert p["kappa"] < 0.05          # near-zero reversion -> very long half-life
    assert p["half_life"] > 20.0


def test_anti_persistent_raises():
    # b < 0 (oscillating) is not a valid OU discretization.
    rng = random.Random(5)
    x = [0.0]
    for _ in range(500):
        x.append(-0.6 * x[-1] + rng.gauss(0, 1))
    with pytest.raises(ValueError):
        fit_ornstein_uhlenbeck(x)


def test_validation():
    with pytest.raises(ValueError):
        fit_ornstein_uhlenbeck([1.0, 2.0])
    with pytest.raises(ValueError):
        fit_ornstein_uhlenbeck([1.0, 2.0, 3.0], dt=0.0)
