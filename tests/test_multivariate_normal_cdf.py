"""Bivariate and trivariate normal CDFs."""

import math

import pytest

from quantforge import bivariate_normal_cdf, trivariate_normal_cdf
from quantforge.mathfns import norm_cdf


def test_bvn_zero_correlation_is_product():
    assert abs(bivariate_normal_cdf(0.5, 0.3, 0.0)
               - norm_cdf(0.5) * norm_cdf(0.3)) < 1e-9


def test_bvn_orthant_formula():
    # P(X<0, Y<0; rho) = 1/4 + asin(rho)/(2 pi).
    for rho in (0.0, 0.5, -0.5, 0.8):
        exact = 0.25 + math.asin(rho) / (2 * math.pi)
        assert abs(bivariate_normal_cdf(0, 0, rho) - exact) < 1e-7


def test_tvn_zero_correlation_is_product():
    got = trivariate_normal_cdf(0.5, 0.3, 0.2, 0, 0, 0)
    assert abs(got - norm_cdf(0.5) * norm_cdf(0.3) * norm_cdf(0.2)) < 1e-9


def test_tvn_reduces_to_bvn_when_third_infinite():
    tv = trivariate_normal_cdf(0.5, 0.3, 8.0, 0.4, 0.2, 0.3)
    bv = bivariate_normal_cdf(0.5, 0.3, 0.4)
    assert abs(tv - bv) < 1e-4


def test_tvn_monotone_in_correlation():
    vals = [trivariate_normal_cdf(0, 0, 0, r, r, r) for r in (0.0, 0.3, 0.6)]
    assert vals[0] < vals[1] < vals[2]
    assert abs(vals[0] - 0.125) < 1e-4        # independent -> 0.5^3


def test_tvn_permutation_symmetry():
    t1 = trivariate_normal_cdf(0.5, 1.0, -0.5, 0.3, 0.4, 0.2)
    t2 = trivariate_normal_cdf(1.0, 0.5, -0.5, 0.3, 0.2, 0.4)
    assert abs(t1 - t2) < 1e-6


def test_bounds():
    v = trivariate_normal_cdf(-1, 0.5, 2, 0.5, -0.3, 0.2)
    assert 0.0 <= v <= 1.0


def test_validation():
    with pytest.raises(ValueError):
        bivariate_normal_cdf(0, 0, 1.5)
    with pytest.raises(ValueError):
        trivariate_normal_cdf(0, 0, 0, 1.5, 0, 0)
