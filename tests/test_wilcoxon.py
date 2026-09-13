"""Wilcoxon signed-rank test and sign test."""

import random

import pytest

from quantforge import wilcoxon_signed_rank_test, sign_test
from quantforge.distributions import binomial_cdf


def test_w_plus_statistic():
    # abs = [1,2,3,4,5], positives at values 1,3,4 -> ranks 1+3+4 = 8.
    assert wilcoxon_signed_rank_test([1, -2, 3, 4, -5])["statistic"] == 8.0


def test_all_positive_textbook():
    d = [1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30]
    r = wilcoxon_signed_rank_test(d)
    assert r["statistic"] == 45.0             # sum of ranks 1..9
    assert abs(r["z"] - 2.6063) < 1e-3
    assert abs(r["p_value"] - 0.00915) < 1e-3


def test_paired_equals_one_sample_on_diffs():
    a = [5, 3, 8, 6, 7]
    b = [4, 4, 6, 5, 5]
    r1 = wilcoxon_signed_rank_test(a, y=b)
    r2 = wilcoxon_signed_rank_test([a[i] - b[i] for i in range(5)])
    assert r1["statistic"] == r2["statistic"]


def test_zeros_dropped():
    assert wilcoxon_signed_rank_test([0, 0, 1, 2, 3])["n"] == 3


def test_symmetric_not_significant_shift_significant():
    rng = random.Random(3)
    sym = [rng.gauss(0, 1) for _ in range(50)]
    assert wilcoxon_signed_rank_test(sym)["p_value"] > 0.05
    shifted = [rng.gauss(1.0, 1) for _ in range(50)]
    assert wilcoxon_signed_rank_test(shifted)["p_value"] < 0.01


def test_sign_test_exact_binomial():
    r = sign_test([1, 1, 1, 1, 1, 1, 1, 1, -1, -1])   # n_plus=8, n=10
    assert r["n_plus"] == 8
    assert abs(r["p_value"] - 2 * binomial_cdf(2, 10, 0.5)) < 1e-12


def test_sign_test_all_positive():
    assert sign_test([1, 2, 3, 4, 5, 6, 7, 8])["p_value"] < 0.01


def test_empty_after_drop_is_pvalue_one():
    assert wilcoxon_signed_rank_test([0, 0, 0])["p_value"] == 1.0
    assert sign_test([0, 0])["p_value"] == 1.0


def test_validation():
    with pytest.raises(ValueError):
        wilcoxon_signed_rank_test([1, 2], y=[1])
    with pytest.raises(ValueError):
        sign_test([1, 2], y=[1])
