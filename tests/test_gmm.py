"""1-D Gaussian mixture model EM."""

import math
import random

import pytest

from quantforge import fit_gaussian_mixture


def test_recovers_separated_components():
    rng = random.Random(1)
    data = [rng.gauss(-3, 0.7) if rng.random() < 0.5 else rng.gauss(3, 0.7)
            for _ in range(2000)]
    r = fit_gaussian_mixture(data, k=2)
    means = sorted(r["means"])
    assert abs(means[0] + 3) < 0.2
    assert abs(means[1] - 3) < 0.2


def test_weights_sum_to_one():
    rng = random.Random(2)
    data = [rng.gauss(0, 1) if rng.random() < 0.3 else rng.gauss(5, 1)
            for _ in range(1500)]
    r = fit_gaussian_mixture(data, k=2)
    assert abs(sum(r["weights"]) - 1.0) < 1e-9
    assert all(w >= 0.0 for w in r["weights"])


def test_single_component_is_sample_moments():
    rng = random.Random(3)
    data = [rng.gauss(5, 2) for _ in range(1000)]
    r = fit_gaussian_mixture(data, k=1)
    m = sum(data) / len(data)
    v = sum((x - m) ** 2 for x in data) / len(data)
    assert abs(r["means"][0] - m) < 1e-6
    assert abs(r["variances"][0] - v) < 1e-6


def test_log_likelihood_non_decreasing():
    rng = random.Random(4)
    data = [rng.gauss(-2, 1) if rng.random() < 0.6 else rng.gauss(2, 1)
            for _ in range(1000)]
    lls = [fit_gaussian_mixture(data, k=2, max_iter=m)["log_likelihood"]
           for m in (1, 2, 3, 5, 10, 50)]
    assert all(lls[i] <= lls[i + 1] + 1e-6 for i in range(len(lls) - 1))


def test_recovers_mixing_weights():
    rng = random.Random(5)
    data = [rng.gauss(-4, 0.6) if rng.random() < 0.7 else rng.gauss(4, 0.6)
            for _ in range(3000)]
    r = fit_gaussian_mixture(data, k=2)
    # The component near -4 should carry ~70% weight.
    lo_idx = 0 if r["means"][0] < r["means"][1] else 1
    assert abs(r["weights"][lo_idx] - 0.7) < 0.05


def test_validation():
    with pytest.raises(ValueError):
        fit_gaussian_mixture([1.0], k=2)
    with pytest.raises(ValueError):
        fit_gaussian_mixture([2.0, 2.0, 2.0], k=1)   # zero spread
