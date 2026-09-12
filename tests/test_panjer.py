"""Panjer recursion for compound-Poisson aggregate loss."""

import math

import pytest

from quantforge import (
    panjer_poisson, aggregate_mean, stop_loss_premium, layer_expected_loss,
)


SEV = [0.0, 0.4, 0.6]      # E[X]=1.6, E[X^2]=2.8
LAM = 3.0


def test_distribution_sums_to_one():
    g = panjer_poisson(LAM, SEV)
    assert abs(sum(g) - 1.0) < 1e-6


def test_mean_and_variance_identities():
    g = panjer_poisson(LAM, SEV)
    assert abs(aggregate_mean(g) - LAM * 1.6) < 1e-4
    var = sum((k - aggregate_mean(g)) ** 2 * g[k] for k in range(len(g)))
    assert abs(var - LAM * 2.8) < 0.01


def test_g0_and_g1_closed_form():
    g = panjer_poisson(LAM, SEV)
    assert abs(g[0] - math.exp(-LAM)) < 1e-9
    assert abs(g[1] - math.exp(-LAM) * LAM * 0.4) < 1e-9


def test_stop_loss_at_zero_is_mean_and_decreasing():
    g = panjer_poisson(LAM, SEV)
    assert abs(stop_loss_premium(g, 0) - aggregate_mean(g)) < 1e-4
    assert stop_loss_premium(g, 2) > stop_loss_premium(g, 5)


def test_layer_full_equals_mean_and_capped():
    g = panjer_poisson(LAM, SEV)
    assert abs(layer_expected_loss(g, 0, 1e9) - aggregate_mean(g)) < 1e-4
    assert layer_expected_loss(g, 3, 2) <= 2.0


def test_validation():
    with pytest.raises(ValueError):
        panjer_poisson(-1, SEV)
    with pytest.raises(ValueError):
        panjer_poisson(LAM, [])
