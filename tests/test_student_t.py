"""Student's t distribution and fat-tailed parametric VaR / ES."""

import pytest

import math
import random

from quantforge import (
    t_pdf, t_cdf, t_ppf, student_t_var, student_t_expected_shortfall,
    fit_df_from_kurtosis, fit_student_t,
)
from quantforge.mathfns import norm_cdf, norm_ppf


def test_cdf_center_and_symmetry():
    assert t_cdf(0, 5) == pytest.approx(0.5)
    assert t_cdf(-1.5, 5) == pytest.approx(1 - t_cdf(1.5, 5), abs=1e-9)


def test_converges_to_normal():
    assert t_cdf(1.0, 1e6) == pytest.approx(norm_cdf(1.0), abs=1e-4)
    assert t_ppf(0.95, 1e6) == pytest.approx(norm_ppf(0.95), abs=1e-3)


def test_critical_value_df5():
    assert t_ppf(0.95, 5) == pytest.approx(2.015, abs=0.01)


def test_ppf_symmetry():
    assert t_ppf(0.05, 10) == pytest.approx(-t_ppf(0.95, 10), abs=1e-6)


def test_pdf_integrates_to_one():
    xs = [x * 0.01 for x in range(-1500, 1501)]
    integ = sum(t_pdf(x, 5) * 0.01 for x in xs)
    assert integ == pytest.approx(1.0, abs=1e-3)


def test_t_var_fatter_than_normal():
    nvar = -(0.0 + 0.02 * norm_ppf(0.01))
    assert student_t_var(0.0, 0.02, 4, 0.99) > nvar


def test_es_at_least_var():
    assert student_t_expected_shortfall(0.0, 0.02, 4, 0.99) > \
        student_t_var(0.0, 0.02, 4, 0.99)


def test_var_converges_to_normal():
    nvar = -(0.0 + 0.02 * norm_ppf(0.01))
    assert student_t_var(0.0, 0.02, 1e6, 0.99) == pytest.approx(nvar, abs=1e-4)


def test_df_from_kurtosis_formula():
    # Excess kurtosis 6/(df-4); df=8 -> ek=1.5.
    assert fit_df_from_kurtosis(1.5) == pytest.approx(8.0)
    assert fit_df_from_kurtosis(3.0) < fit_df_from_kurtosis(1.0)   # higher kurt -> lower df


def test_fit_recovers_df_order():
    random.seed(7)

    def t_sample(df):
        z = random.gauss(0, 1)
        chi = sum(random.gauss(0, 1) ** 2 for _ in range(df))
        return z / math.sqrt(chi / df)

    data = [t_sample(8) for _ in range(50000)]
    _, _, df = fit_student_t(data)
    assert 4 < df < 20   # noisy moment match, but finite and near 8


def test_fit_validation():
    with pytest.raises(ValueError):
        fit_df_from_kurtosis(-1)
    with pytest.raises(ValueError):
        fit_student_t([1, 2, 3])


def test_validation():
    with pytest.raises(ValueError):
        t_pdf(0, 0)
    with pytest.raises(ValueError):
        student_t_expected_shortfall(0.0, 0.02, 1, 0.99)   # df <= 1
