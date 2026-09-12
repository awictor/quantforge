"""Buhlmann and Buhlmann-Straub credibility."""

import pytest

from quantforge import (
    buhlmann_k, credibility_factor, buhlmann_premium, buhlmann_straub_premium,
)


def test_k_and_z_bounds():
    k = buhlmann_k(100, 20)
    assert k == 5.0
    assert 0.0 <= credibility_factor(10, k) <= 1.0


def test_premium_between_own_and_collective():
    p = buhlmann_premium(500, 400, 10, 100, 20)
    assert 400 < p < 500


def test_credibility_limits():
    k = buhlmann_k(100, 20)
    assert credibility_factor(1e9, k) > 0.999      # n -> inf
    assert credibility_factor(0, k) == 0.0          # n = 0
    assert abs(buhlmann_premium(500, 400, 0, 100, 20) - 400) < 1e-9


def test_variance_ratio_effects():
    assert credibility_factor(10, buhlmann_k(200, 20)) < credibility_factor(10, buhlmann_k(100, 20))
    assert credibility_factor(10, buhlmann_k(100, 40)) > credibility_factor(10, buhlmann_k(100, 20))


def test_straub_matches_buhlmann_on_equal_exposure():
    bs = buhlmann_straub_premium([50, 50, 50], [10, 10, 10], 4.0, 100, 20)
    assert abs(bs - buhlmann_premium(5.0, 4.0, 30, 100, 20)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        buhlmann_k(100, 0)
    with pytest.raises(ValueError):
        buhlmann_straub_premium([1], [1, 2], 4, 100, 20)
