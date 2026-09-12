"""Student-t copula sampler."""

import pytest

from quantforge import student_t_copula_sample, gaussian_copula_sample
from quantforge.copula_stats import spearman_rho
from quantforge.tail_dependence import lower_tail_dependence


def test_rank_correlation_near_target():
    s = student_t_copula_sample([[1.0, 0.6], [0.6, 1.0]], df=5, n=5000, seed=42)
    x = [r[0] for r in s]
    y = [r[1] for r in s]
    assert abs(spearman_rho(x, y) - 0.6) < 0.08


def test_margins_uniform():
    s = student_t_copula_sample([[1.0, 0.6], [0.6, 1.0]], df=5, n=5000, seed=42)
    x = [r[0] for r in s]
    assert abs(sum(x) / len(x) - 0.5) < 0.02
    assert all(0.0 < v < 1.0 for v in x)


def test_lower_df_has_more_tail_dependence():
    R = [[1.0, 0.6], [0.6, 1.0]]
    lo = student_t_copula_sample(R, df=3, n=8000, seed=1)
    hi = student_t_copula_sample(R, df=100, n=8000, seed=1)
    td_lo = lower_tail_dependence([r[0] for r in lo], [r[1] for r in lo], 0.05)
    td_hi = lower_tail_dependence([r[0] for r in hi], [r[1] for r in hi], 0.05)
    assert td_lo > td_hi


def test_high_df_approaches_gaussian_copula():
    R = [[1.0, 0.6], [0.6, 1.0]]
    hi = student_t_copula_sample(R, df=100, n=8000, seed=1)
    g = gaussian_copula_sample(R, 8000, seed=1)
    td_hi = lower_tail_dependence([r[0] for r in hi], [r[1] for r in hi], 0.05)
    td_g = lower_tail_dependence([r[0] for r in g], [r[1] for r in g], 0.05)
    assert abs(td_hi - td_g) < 0.15


def test_independence_near_zero_rank_correlation():
    s = student_t_copula_sample([[1.0, 0.0], [0.0, 1.0]], df=5, n=5000, seed=7)
    assert abs(spearman_rho([r[0] for r in s], [r[1] for r in s])) < 0.05


def test_validation():
    R = [[1.0, 0.6], [0.6, 1.0]]
    with pytest.raises(ValueError):
        student_t_copula_sample(R, df=0, n=10)
    with pytest.raises(ValueError):
        student_t_copula_sample(R, df=5, n=0)
    with pytest.raises(ValueError):
        student_t_copula_sample([[1.0, 1.5], [1.5, 1.0]], df=5, n=10)  # non-PD
