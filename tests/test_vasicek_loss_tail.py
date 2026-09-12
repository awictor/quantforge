"""Vasicek large-pool loss density and expected shortfall."""

import math

import pytest

from quantforge import (vasicek_loss_pdf, vasicek_loss_cdf, vasicek_loss_quantile,
                        vasicek_loss_expected_shortfall)


def test_density_integrates_to_one_and_mean_is_pd():
    pd, rho = 0.02, 0.15
    n = 200000
    dx = 1.0 / n
    total = 0.0
    mean = 0.0
    for k in range(1, n):
        x = k * dx
        f = vasicek_loss_pdf(x, pd, rho)
        total += f * dx
        mean += x * f * dx
    assert abs(total - 1.0) < 1e-3
    assert abs(mean - pd) < 1e-4


def test_density_matches_cdf_finite_difference():
    pd, rho = 0.03, 0.2
    h = 1e-6
    for x in (0.01, 0.05, 0.1):
        fd = (vasicek_loss_cdf(x + h, pd, rho) - vasicek_loss_cdf(x - h, pd, rho)) / (2 * h)
        assert abs(vasicek_loss_pdf(x, pd, rho) - fd) < 1e-4


def test_es_matches_tail_average_of_quantile():
    pd, rho, q = 0.02, 0.15, 0.99
    es = vasicek_loss_expected_shortfall(q, pd, rho)
    m = 200000
    s = sum(vasicek_loss_quantile(q + (1 - q) * (k + 0.5) / m, pd, rho) for k in range(m))
    s /= m
    assert abs(es - s) < 1e-4


def test_es_at_least_var():
    pd, rho = 0.02, 0.15
    for q in (0.95, 0.99, 0.999):
        assert vasicek_loss_expected_shortfall(q, pd, rho) >= vasicek_loss_quantile(q, pd, rho)


def test_es_between_mean_and_one():
    pd, rho, q = 0.05, 0.25, 0.99
    es = vasicek_loss_expected_shortfall(q, pd, rho)
    assert pd < es < 1.0


def test_es_increasing_in_q_and_rho():
    pd = 0.03
    e1 = vasicek_loss_expected_shortfall(0.95, pd, 0.2)
    e2 = vasicek_loss_expected_shortfall(0.99, pd, 0.2)
    assert e2 > e1
    e3 = vasicek_loss_expected_shortfall(0.99, pd, 0.3)
    assert e3 > e2


def test_es_matches_monte_carlo():
    from quantforge.mathfns import norm_ppf
    import random
    pd, rho, q = 0.02, 0.15, 0.99
    thr = norm_ppf(pd)
    rng = random.Random(12345)
    n = 400000
    losses = []
    for _ in range(n):
        m = rng.gauss(0, 1)
        p = 0.5 * (1 + math.erf(((thr - math.sqrt(rho) * m) / math.sqrt(1 - rho)) / math.sqrt(2)))
        losses.append(p)
    losses.sort()
    idx = int(q * n)
    es_mc = sum(losses[idx:]) / (n - idx)
    assert abs(vasicek_loss_expected_shortfall(q, pd, rho) - es_mc) < 5e-4


def test_validation():
    with pytest.raises(ValueError):
        vasicek_loss_pdf(0.0, 0.02, 0.15)
    with pytest.raises(ValueError):
        vasicek_loss_pdf(0.5, 1.0, 0.15)
    with pytest.raises(ValueError):
        vasicek_loss_expected_shortfall(1.0, 0.02, 0.15)
    with pytest.raises(ValueError):
        vasicek_loss_expected_shortfall(0.99, 0.02, 0.0)
