"""Student's t distribution and fat-tailed parametric VaR / ES."""

import pytest

from quantforge import (
    t_pdf, t_cdf, t_ppf, student_t_var, student_t_expected_shortfall,
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


def test_validation():
    with pytest.raises(ValueError):
        t_pdf(0, 0)
    with pytest.raises(ValueError):
        student_t_expected_shortfall(0.0, 0.02, 1, 0.99)   # df <= 1
