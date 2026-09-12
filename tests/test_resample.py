"""Bootstrap and jackknife resampling."""

import statistics

import pytest

from quantforge import bootstrap_ci, stationary_bootstrap_ci, jackknife_estimate


DATA = [3, 5, 7, 4, 6, 5, 8, 2, 5, 6, 4, 7, 5, 3, 6, 5, 4, 7, 6, 5]


def test_bootstrap_ci_brackets_point():
    lo, pt, hi = bootstrap_ci(DATA, n_boot=1000)
    assert lo <= pt <= hi


def test_bootstrap_deterministic():
    assert bootstrap_ci(DATA, n_boot=500) == bootstrap_ci(DATA, n_boot=500)


def test_jackknife_se_matches_analytic():
    est, se = jackknife_estimate(DATA)
    approx = statistics.stdev(DATA) / len(DATA) ** 0.5
    assert se == pytest.approx(approx, abs=0.05)
    assert est == pytest.approx(statistics.mean(DATA))


def test_custom_statistic():
    assert bootstrap_ci(DATA, statistic=max, n_boot=500)[1] == max(DATA)


def test_stationary_bootstrap_deterministic():
    a = stationary_bootstrap_ci(DATA, mean_block=5, n_boot=400)
    b = stationary_bootstrap_ci(DATA, mean_block=5, n_boot=400)
    assert a == b


def test_stationary_ci_brackets_point():
    lo, pt, hi = stationary_bootstrap_ci(DATA, mean_block=5, n_boot=600)
    assert lo <= pt <= hi


def test_validation():
    with pytest.raises(ValueError):
        bootstrap_ci([], n_boot=100)
    with pytest.raises(ValueError):
        jackknife_estimate([5])
    with pytest.raises(ValueError):
        stationary_bootstrap_ci(DATA, mean_block=0)
