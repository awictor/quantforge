"""Exponential-kernel Hawkes self-exciting process."""

import math

import pytest

from quantforge import (hawkes_intensity, hawkes_branching_ratio,
                        hawkes_log_likelihood, hawkes_simulate, hawkes_fit)


def _brute_ll(events, mu, al, be, T):
    ll = 0.0
    for i, ti in enumerate(events):
        a = sum(math.exp(-be * (ti - tj)) for tj in events[:i])
        ll += math.log(mu + al * a)
    comp = mu * T + sum((al / be) * (1 - math.exp(-be * (T - ti))) for ti in events)
    return ll - comp


def test_log_likelihood_recursion_matches_brute():
    ev = [0.5, 1.2, 1.5, 3.0, 3.1, 5.0]
    for mu, al, be in [(0.4, 0.6, 1.5), (1.0, 0.2, 0.9), (0.3, 0.0, 2.0)]:
        assert abs(hawkes_log_likelihood(ev, mu, al, be, 6.0)
                   - _brute_ll(ev, mu, al, be, 6.0)) < 1e-10


def test_branching_ratio():
    assert abs(hawkes_branching_ratio(0.6, 1.5) - 0.4) < 1e-12


def test_intensity_adds_kernel():
    # Intensity after two events = mu + alpha(e^{-beta d1} + e^{-beta d2}).
    hist = [1.0, 2.0]
    mu, al, be = 0.5, 0.8, 1.2
    expected = mu + al * (math.exp(-be * (3.0 - 1.0)) + math.exp(-be * (3.0 - 2.0)))
    assert abs(hawkes_intensity(3.0, hist, mu, al, be) - expected) < 1e-12


def test_poisson_limit_rate():
    # alpha = 0 reduces to a Poisson process with rate mu.
    ev = hawkes_simulate(0.5, 0.0, 2.0, 4000.0, seed=7)
    assert abs(len(ev) / 4000.0 - 0.5) < 0.06


@pytest.mark.slow
def test_stationary_rate_matches_theory():
    # Mean intensity of a stationary Hawkes is mu / (1 - alpha/beta).
    for al, exp_rate in [(0.2, 0.5 / 0.9), (0.5, 0.5 / 0.75), (0.9, 0.5 / 0.55)]:
        ev = hawkes_simulate(0.5, al, 2.0, 10000.0, seed=7)
        assert abs(len(ev) / 10000.0 - exp_rate) < 0.06 * exp_rate


@pytest.mark.slow
def test_mle_recovers_parameters():
    true = (0.5, 0.9, 2.0)
    ev = hawkes_simulate(*true, t_max=10000.0, seed=99)
    r = hawkes_fit(ev)
    assert abs(r["mu"] - 0.5) < 0.1
    assert abs(r["alpha"] - 0.9) < 0.2
    assert abs(r["beta"] - 2.0) < 0.4
    assert 0.0 < r["branching_ratio"] < 1.0


def test_validation():
    with pytest.raises(ValueError):
        hawkes_intensity(1.0, [], 0.0, 0.5, 1.0)
    with pytest.raises(ValueError):
        hawkes_log_likelihood([], 0.5, 0.5, 1.0)
    with pytest.raises(ValueError):
        hawkes_fit([1.0])
