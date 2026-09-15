import math

import pytest

from quantforge import (
    beta_binomial_posterior,
    gamma_poisson_posterior,
    normal_normal_posterior,
    hpd_interval,
)


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _integ(f, a, b, N=20000):
    h = (b - a) / N
    s = 0.5 * (f(a) + f(b))
    for i in range(1, N):
        s += f(a + i * h)
    return s * h


def test_beta_binomial():
    p = beta_binomial_posterior(2, 2, 7, 10)
    assert p["alpha"] == 9 and p["beta"] == 5
    assert close(p["mean"], 9 / 14)
    assert close(p["var"], 9 * 5 / ((14 ** 2) * 15))


def test_laplace_rule():
    p = beta_binomial_posterior(1, 1, 3, 10)
    assert close(p["mean"], 4 / 12)


def test_gamma_poisson():
    p = gamma_poisson_posterior(2, 1, 15, 5)
    assert p["shape"] == 17 and p["rate"] == 6
    assert close(p["mean"], 17 / 6)
    assert close(p["var"], 17 / 36)


def test_normal_normal_precision_weighting():
    p = normal_normal_posterior(0.0, 1.0, [2.0, 2.0, 2.0, 2.0], 1.0)
    assert close(p["var"], 0.2)
    assert close(p["mean"], 1.6)


def test_vague_prior_recovers_mle():
    p = normal_normal_posterior(0.0, 1e12, [2.0, 2.0, 2.0, 2.0], 1.0)
    assert close(p["mean"], 2.0, 1e-6)


def test_hpd_normal_matches_analytic():
    mu, sig = 1.6, math.sqrt(0.2)
    pdf = lambda x: math.exp(-0.5 * ((x - mu) / sig) ** 2)
    lo, hi = hpd_interval(pdf, mu - 8 * sig, mu + 8 * sig, mass=0.95, n_grid=20000)
    assert close(lo, mu - 1.96 * sig, 2e-2)
    assert close(hi, mu + 1.96 * sig, 2e-2)
    assert close((lo + hi) / 2, mu, 1e-2)


def test_hpd_beta_contains_mean_and_correct_mass():
    a, b = 9, 5
    pdf = lambda x: x ** (a - 1) * (1 - x) ** (b - 1) if 0 < x < 1 else 0.0
    lo, hi = hpd_interval(pdf, 0.0, 1.0, mass=0.95, n_grid=20000)
    assert lo < a / (a + b) < hi
    tot = _integ(pdf, 0.0, 1.0)
    inside = _integ(pdf, lo, hi)
    assert abs(inside / tot - 0.95) < 0.03


def test_errors():
    with pytest.raises(ValueError):
        beta_binomial_posterior(0, 2, 1, 10)
    with pytest.raises(ValueError):
        beta_binomial_posterior(2, 2, 11, 10)
    with pytest.raises(ValueError):
        gamma_poisson_posterior(-1, 1, 5, 5)
    with pytest.raises(ValueError):
        normal_normal_posterior(0.0, -1.0, [1.0], 1.0)
    with pytest.raises(ValueError):
        hpd_interval(lambda x: 1.0, 0.0, 1.0, mass=1.5)
