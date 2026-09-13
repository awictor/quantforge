"""Statistical power and sample-size calculations."""

import math
import random

import pytest

from quantforge import (two_sample_t_power, two_sample_t_sample_size,
                        one_sample_z_power, one_sample_z_sample_size,
                        proportion_power, proportion_sample_size)
from quantforge.hypothesis import two_sample_t_test


def test_two_sample_sample_size_textbook():
    # d = 0.5, power 0.8, alpha 0.05 -> ~63-64 per group.
    n = two_sample_t_sample_size(0.5, 0.80)
    assert 60 <= n <= 66


def test_sample_size_delivers_target_power():
    for d in (0.3, 0.5, 0.8):
        n = two_sample_t_sample_size(d, 0.80)
        assert two_sample_t_power(d, n) >= 0.80


def test_power_matches_monte_carlo():
    d, n, alpha = 0.6, 40, 0.05
    computed = two_sample_t_power(d, n, alpha)
    rng = random.Random(7)
    rej, trials = 0, 4000
    for _ in range(trials):
        a = [rng.gauss(d, 1) for _ in range(n)]
        b = [rng.gauss(0, 1) for _ in range(n)]
        _, p = two_sample_t_test(a, b, equal_var=True)
        if p < alpha:
            rej += 1
    assert abs(computed - rej / trials) < 0.03


def test_power_monotone_in_effect_and_n():
    assert two_sample_t_power(0.3, 40) < two_sample_t_power(0.8, 40)
    assert two_sample_t_power(0.5, 20) < two_sample_t_power(0.5, 100)


def test_zero_effect_power_equals_alpha():
    assert abs(two_sample_t_power(1e-12, 40, 0.05) - 0.05) < 1e-3


def test_one_sample_sample_size_textbook():
    n = one_sample_z_sample_size(0.5, 0.80)
    assert 30 <= n <= 34
    assert one_sample_z_power(0.5, n) >= 0.80


def test_proportion_sample_size_textbook():
    # p1=0.5, p2=0.65, power 0.8 -> ~169-170 per group.
    n = proportion_sample_size(0.5, 0.65, 0.80)
    assert 165 <= n <= 175
    assert proportion_power(0.5, 0.65, n) >= 0.80


def test_proportion_power_monotone():
    assert proportion_power(0.5, 0.55, 100) < proportion_power(0.5, 0.70, 100)


def test_validation():
    with pytest.raises(ValueError):
        two_sample_t_power(0.5, 1)
    with pytest.raises(ValueError):
        two_sample_t_sample_size(0.0)
    with pytest.raises(ValueError):
        two_sample_t_sample_size(0.5, power=1.0)
    with pytest.raises(ValueError):
        proportion_power(0.0, 0.5, 50)
    with pytest.raises(ValueError):
        proportion_sample_size(0.5, 0.5)
