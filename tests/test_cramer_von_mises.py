"""Two-sample Cramer-von Mises test."""

import random

import pytest

from quantforge import cramer_von_mises_2samp
from quantforge.cramer_von_mises import _cvm_asymptotic_sf


def test_asymptotic_tail_matches_critical_values():
    # Published upper-tail critical values of the limiting CvM distribution.
    assert abs(_cvm_asymptotic_sf(0.461) - 0.05) < 0.005
    assert abs(_cvm_asymptotic_sf(0.743) - 0.01) < 0.005
    assert abs(_cvm_asymptotic_sf(0.347) - 0.10) < 0.005


def test_identical_samples_zero_statistic():
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert abs(cramer_von_mises_2samp(a, a)["statistic"]) < 1e-9


def test_agrees_with_permutation_pvalue():
    def perm_p(a, b, nperm=2000, seed=1):
        rng = random.Random(seed)
        obs = cramer_von_mises_2samp(a, b)["statistic"]
        pool = list(a) + list(b)
        n = len(a)
        cnt = 0
        for _ in range(nperm):
            rng.shuffle(pool)
            s = cramer_von_mises_2samp(pool[:n], pool[n:])["statistic"]
            if s >= obs - 1e-12:
                cnt += 1
        return (1 + cnt) / (1 + nperm)

    rng = random.Random(5)
    a = [rng.gauss(0, 1) for _ in range(40)]
    b = [rng.gauss(0, 1) for _ in range(40)]
    r = cramer_von_mises_2samp(a, b)
    assert abs(r["p_value"] - perm_p(a, b)) < 0.06


def test_detects_shift():
    rng = random.Random(5)
    a = [rng.gauss(0, 1) for _ in range(40)]
    b = [rng.gauss(1.5, 1) for _ in range(40)]
    assert cramer_von_mises_2samp(a, b)["p_value"] < 0.01


def test_null_rejection_rate_near_alpha():
    rng = random.Random(9)
    rej = 0
    N = 400
    for _ in range(N):
        a = [rng.gauss(0, 1) for _ in range(30)]
        b = [rng.gauss(0, 1) for _ in range(30)]
        if cramer_von_mises_2samp(a, b)["p_value"] < 0.05:
            rej += 1
    assert 0.01 <= rej / N <= 0.09


def test_validation():
    with pytest.raises(ValueError):
        cramer_von_mises_2samp([], [1, 2])
    with pytest.raises(ValueError):
        cramer_von_mises_2samp([1, 2], [])
