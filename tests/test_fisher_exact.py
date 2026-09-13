"""Fisher's exact test for 2x2 tables."""

from math import comb

import pytest

from quantforge import fisher_exact_test


def test_tea_tasting_two_sided():
    r = fisher_exact_test([[3, 1], [1, 3]])
    assert abs(r["odds_ratio"] - 9.0) < 1e-9
    assert abs(r["p_value"] - 0.4857) < 1e-3


def test_tea_tasting_one_sided():
    assert abs(fisher_exact_test([[3, 1], [1, 3]], "greater")["p_value"] - 0.2429) < 1e-3


def test_known_table_matches_scipy():
    # scipy.stats.fisher_exact([[1,9],[11,3]]): OR=0.0303, one-sided less p=0.00135.
    r = fisher_exact_test([[1, 9], [11, 3]])
    assert abs(r["odds_ratio"] - 0.030303) < 1e-4
    less = fisher_exact_test([[1, 9], [11, 3]], "less")["p_value"]
    assert abs(less - 0.00135) < 5e-4


def test_probabilities_sum_to_one():
    a, b, c, d = 3, 1, 1, 3
    r1, r2, c1, n = a + b, c + d, a + c, a + b + c + d
    lo, hi = max(0, c1 - r2), min(r1, c1)
    total = sum(comb(r1, k) * comb(r2, c1 - k) / comb(n, c1) for k in range(lo, hi + 1))
    assert abs(total - 1.0) < 1e-9


def test_strong_association_small_p():
    assert fisher_exact_test([[20, 2], [3, 25]])["p_value"] < 1e-6


def test_independence_p_near_one():
    r = fisher_exact_test([[10, 10], [10, 10]])
    assert abs(r["p_value"] - 1.0) < 1e-9
    assert abs(r["odds_ratio"] - 1.0) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        fisher_exact_test([[1, 2, 3], [4, 5, 6]])
    with pytest.raises(ValueError):
        fisher_exact_test([[1, -1], [2, 3]])
    with pytest.raises(ValueError):
        fisher_exact_test([[1, 2], [3, 4]], alternative="bad")
