"""Negative-binomial Panjer recursion (over-dispersed aggregate loss)."""

import pytest

from quantforge import panjer_negative_binomial as nb, panjer_poisson as pp, aggregate_mean as am


SEV = [0.0, 0.4, 0.6]      # E[X]=1.6, E[X^2]=2.8


def test_sums_to_one():
    g = nb(5, 0.5, SEV)
    assert abs(sum(g) - 1.0) < 1e-6


def test_mean_and_variance_identities():
    g = nb(5, 0.5, SEV)                 # mean N = 5, var N = 10
    assert abs(am(g) - 5 * 1.6) < 1e-3
    var = sum((k - am(g)) ** 2 * g[k] for k in range(len(g)))
    # E[N]Var[X] + Var[N]E[X]^2 = 5*0.24 + 10*2.56 = 26.8
    assert abs(var - 26.8) < 0.1


def test_over_dispersed_versus_poisson():
    g_nb = nb(5, 0.5, SEV)              # mean S = 8
    g_p = pp(5.0, SEV)                  # mean S = 8, var = 14
    var_nb = sum((k - am(g_nb)) ** 2 * g_nb[k] for k in range(len(g_nb)))
    var_p = sum((k - am(g_p)) ** 2 * g_p[k] for k in range(len(g_p)))
    assert abs(am(g_nb) - am(g_p)) < 1e-3      # same mean
    assert var_nb > var_p                       # heavier tail


def test_validation():
    with pytest.raises(ValueError):
        nb(-1, 0.5, SEV)
    with pytest.raises(ValueError):
        nb(5, 1.5, SEV)
    with pytest.raises(ValueError):
        nb(5, 0.5, [])
