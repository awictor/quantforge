"""Divergences between discrete distributions."""

import math
import random

import pytest

from quantforge import (kl_divergence, jensen_shannon_divergence, hellinger_distance,
                        total_variation_distance, bhattacharyya_distance)


def test_identical_all_zero():
    p = [0.2, 0.3, 0.5]
    assert kl_divergence(p, p) == 0.0
    assert jensen_shannon_divergence(p, p) == 0.0
    assert hellinger_distance(p, p) == 0.0
    assert total_variation_distance(p, p) == 0.0
    assert abs(bhattacharyya_distance(p, p)) < 1e-12


def test_kl_asymmetric_js_symmetric():
    p, q = [0.7, 0.3], [0.4, 0.6]
    assert abs(kl_divergence(p, q) - kl_divergence(q, p)) > 1e-6
    assert abs(jensen_shannon_divergence(p, q) - jensen_shannon_divergence(q, p)) < 1e-15


def test_js_bounded_by_log2():
    rng = random.Random(3)
    for _ in range(1000):
        n = rng.randint(2, 6)
        a = [rng.random() for _ in range(n)]
        b = [rng.random() for _ in range(n)]
        assert 0 <= jensen_shannon_divergence(a, b) <= math.log(2) + 1e-12


def test_hellinger_tv_in_unit_interval():
    rng = random.Random(5)
    for _ in range(1000):
        n = rng.randint(2, 6)
        a = [rng.random() for _ in range(n)]
        b = [rng.random() for _ in range(n)]
        assert 0 <= hellinger_distance(a, b) <= 1 + 1e-12
        assert 0 <= total_variation_distance(a, b) <= 1 + 1e-12


def test_disjoint_supports():
    p, q = [1.0, 0.0], [0.0, 1.0]
    assert total_variation_distance(p, q) == 1.0
    assert abs(hellinger_distance(p, q) - 1.0) < 1e-12
    assert bhattacharyya_distance(p, q) == float("inf")


def test_kl_raises_on_zero_support():
    with pytest.raises(ValueError):
        kl_divergence([0.5, 0.5], [1.0, 0.0])


def test_pinsker_inequality():
    rng = random.Random(7)
    for _ in range(200):
        n = rng.randint(2, 5)
        a = [rng.random() + 0.1 for _ in range(n)]
        b = [rng.random() + 0.1 for _ in range(n)]
        tv = total_variation_distance(a, b)
        kl = kl_divergence(a, b)
        assert tv <= math.sqrt(kl / 2) + 1e-9


def test_normalization():
    assert abs(kl_divergence([2, 3, 5], [1, 1, 1])
               - kl_divergence([0.2, 0.3, 0.5], [1 / 3, 1 / 3, 1 / 3])) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        kl_divergence([0.5], [0.5, 0.5])
    with pytest.raises(ValueError):
        hellinger_distance([1.0, -1.0], [0.5, 0.5])
