"""Hawkes time-rescaling residuals and goodness-of-fit."""

import math

import pytest

from quantforge import (hawkes_residuals, hawkes_gof_test, hawkes_simulate,
                        hawkes_fit)


def _brute_tau(events, mu, al, be):
    taus = []
    for i in range(1, len(events)):
        a, b = events[i - 1], events[i]
        comp = mu * (b - a)
        for tj in events:
            if tj <= a:
                comp += (al / be) * (math.exp(-be * (a - tj)) - math.exp(-be * (b - tj)))
        taus.append(comp)
    return taus


def test_residuals_match_brute_compensator():
    ev = [0.5, 1.2, 1.5, 3.0, 3.1, 5.0]
    mu, al, be = 0.4, 0.6, 1.5
    r = hawkes_residuals(ev, mu, al, be)
    bt = _brute_tau(ev, mu, al, be)
    assert max(abs(r[i] - bt[i]) for i in range(len(r))) < 1e-12


def test_residual_count():
    ev = [0.1, 0.5, 0.9, 2.0]
    assert len(hawkes_residuals(ev, 0.5, 0.3, 1.0)) == len(ev) - 1


@pytest.mark.slow
def test_correct_model_residuals_are_unit_exponential():
    ev = hawkes_simulate(0.5, 0.9, 2.0, 6000.0, seed=42)
    f = hawkes_fit(ev)
    res = hawkes_residuals(ev, f["mu"], f["alpha"], f["beta"])
    mean = sum(res) / len(res)
    assert abs(mean - 1.0) < 0.05          # Exp(1) has mean 1


@pytest.mark.slow
def test_gof_accepts_true_model_rejects_wrong():
    ev = hawkes_simulate(0.5, 0.9, 2.0, 6000.0, seed=42)
    f = hawkes_fit(ev)
    _, p_true = hawkes_gof_test(ev, f["mu"], f["alpha"], f["beta"])
    _, p_wrong = hawkes_gof_test(ev, 0.5, 0.1, 5.0)
    assert p_true > 0.05
    assert p_wrong < 0.01


def test_validation():
    with pytest.raises(ValueError):
        hawkes_residuals([1.0], 0.5, 0.3, 1.0)
    with pytest.raises(ValueError):
        hawkes_gof_test([1.0, 2.0], 0.5, 0.3, 1.0)   # < 3 events -> < 2 residuals
