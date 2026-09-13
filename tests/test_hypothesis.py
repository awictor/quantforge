"""Classical hypothesis tests on the distribution CDFs."""

import pytest

from quantforge import (chi_square_gof_test, chi_square_independence_test,
                        one_way_anova, two_sample_t_test, binomial_test,
                        one_sample_t_test, paired_t_test, mann_whitney_u)


def test_gof_uniform_die():
    obs = [16, 18, 16, 14, 12, 12]
    stat, p = chi_square_gof_test(obs)
    assert abs(stat - 2.0) < 1e-9           # by hand: sum (O-14.667)^2/14.667
    assert 0.8 < p < 0.9                     # 5 df, small statistic -> large p


def test_gof_rejects_bad_fit():
    obs = [50, 10, 10, 10, 10, 10]
    stat, p = chi_square_gof_test(obs)
    assert p < 0.001


def test_gof_with_explicit_expected():
    stat, p = chi_square_gof_test([20, 30], [25, 25])
    assert abs(stat - 2.0) < 1e-9           # (25/25 + 25/25) = 2
    assert 0 < p < 1


def test_independence_2x2():
    stat, p = chi_square_independence_test([[10, 20], [30, 40]])
    assert abs(stat - 0.793651) < 1e-4
    assert 0.3 < p < 0.45


def test_independence_detects_dependence():
    stat, p = chi_square_independence_test([[90, 10], [10, 90]])
    assert p < 1e-10


def test_anova_equals_t_squared_for_two_groups():
    a = [5.1, 4.9, 5.5, 5.0, 5.2]
    b = [6.1, 5.9, 6.3, 6.0, 5.8]
    f, pf = one_way_anova(a, b)
    t, pt = two_sample_t_test(a, b, equal_var=True)
    assert abs(f - t * t) < 1e-9
    assert abs(pf - pt) < 1e-9


def test_anova_rejects_unequal_means():
    f, p = one_way_anova([1, 2, 3], [10, 11, 12], [20, 21, 22])
    assert p < 1e-4


def test_anova_equal_means_high_p():
    f, p = one_way_anova([1, 2, 3, 2], [2, 3, 1, 2], [3, 1, 2, 2])
    assert p > 0.5


def test_welch_t_test():
    a = [27, 25, 30, 28, 26]
    b = [20, 22, 19, 21, 23, 18]
    t, p = two_sample_t_test(a, b, equal_var=False)
    assert t > 0
    assert p < 0.001


def test_pooled_and_welch_agree_for_equal_sizes_and_var():
    a = [10, 12, 11, 13, 9]
    b = [14, 16, 15, 17, 13]
    t1, _ = two_sample_t_test(a, b, equal_var=True)
    t2, _ = two_sample_t_test(a, b, equal_var=False)
    assert abs(t1 - t2) < 1e-9      # equal n and similar variance -> same t


def test_binomial_two_sided_symmetric():
    prop, p = binomial_test(8, 10, 0.5)
    assert prop == 0.8
    assert abs(p - 0.109375) < 1e-6


def test_binomial_one_sided():
    _, pg = binomial_test(8, 10, 0.5, "greater")
    _, pl = binomial_test(8, 10, 0.5, "less")
    assert abs(pg - 0.0546875) < 1e-7
    assert pl > 0.9


def test_binomial_certain_when_all_success():
    _, p = binomial_test(10, 10, 0.5, "greater")
    assert abs(p - 0.5 ** 10) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        chi_square_gof_test([5])
    with pytest.raises(ValueError):
        chi_square_independence_test([[1, 2, 3]])
    with pytest.raises(ValueError):
        one_way_anova([1, 2, 3])
    with pytest.raises(ValueError):
        two_sample_t_test([1.0], [2.0, 3.0])
    with pytest.raises(ValueError):
        binomial_test(12, 10, 0.5)
    with pytest.raises(ValueError):
        binomial_test(5, 10, 0.5, "sideways")


def test_one_sample_t_matches_manual():
    import statistics
    x = [5.1, 5.3, 4.9, 5.4, 5.2, 5.0, 5.5]
    t, p = one_sample_t_test(x, 5.0)
    m, sd = statistics.mean(x), statistics.stdev(x)
    assert abs(t - (m - 5.0) / (sd / len(x) ** 0.5)) < 1e-9
    assert 0.0 < p < 1.0


def test_one_sample_t_zero_at_mean():
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    t, p = one_sample_t_test(x, 3.0)   # mu0 = sample mean
    assert abs(t) < 1e-12
    assert abs(p - 1.0) < 1e-9


def test_paired_equals_one_sample_on_diffs():
    a = [10, 12, 14, 11, 13]
    b = [9, 11, 12, 10, 13]
    tp, pp = paired_t_test(a, b)
    td, pd = one_sample_t_test([ai - bi for ai, bi in zip(a, b)], 0.0)
    assert abs(tp - td) < 1e-12
    assert abs(pp - pd) < 1e-12


def test_mann_whitney_disjoint_groups():
    u, p = mann_whitney_u([1, 2, 3, 4], [5, 6, 7, 8])
    assert u == 0.0            # complete separation
    assert p < 0.05


def test_mann_whitney_u_symmetry():
    # U_a + U_b = n_a * n_b, and the reported u is the smaller of the two.
    a = [19, 22, 16, 29, 24]
    b = [20, 11, 17, 12]
    u, p = mann_whitney_u(a, b)
    assert u <= len(a) * len(b) / 2.0
    assert abs(p - 0.1113) < 0.02


def test_mann_whitney_identical_not_significant():
    u, p = mann_whitney_u([1, 2, 3, 4, 5], [1, 2, 3, 4, 5])
    assert abs(p - 1.0) < 1e-9


def test_extra_validation():
    with pytest.raises(ValueError):
        one_sample_t_test([1.0])
    with pytest.raises(ValueError):
        one_sample_t_test([2.0, 2.0, 2.0])       # zero variance
    with pytest.raises(ValueError):
        paired_t_test([1.0, 2.0], [1.0])         # length mismatch
    with pytest.raises(ValueError):
        mann_whitney_u([], [1.0])
