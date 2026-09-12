"""VaR / TVaR on a discrete aggregate loss distribution."""

import pytest

from quantforge import (
    panjer_poisson, aggregate_var, aggregate_tvar, aggregate_mean,
)


def _g():
    return panjer_poisson(3.0, [0.0, 0.4, 0.6])


def _cdf(g, k):
    return sum(g[j] for j in range(k + 1))


def test_var_is_the_quantile():
    g = _g()
    v = aggregate_var(g, 0.95)
    assert _cdf(g, v) >= 0.95
    assert _cdf(g, v - 1) < 0.95


def test_tvar_at_least_var():
    g = _g()
    assert aggregate_tvar(g, 0.95) >= aggregate_var(g, 0.95)


def test_both_monotone_in_confidence():
    g = _g()
    assert aggregate_var(g, 0.99) >= aggregate_var(g, 0.95) >= aggregate_var(g, 0.90)
    assert aggregate_tvar(g, 0.99) >= aggregate_tvar(g, 0.95)


def test_tvar_low_confidence_is_mean():
    g = _g()
    assert abs(aggregate_tvar(g, 0.001) - aggregate_mean(g)) < 0.5


def test_validation():
    with pytest.raises(ValueError):
        aggregate_var(_g(), 1.0)
    with pytest.raises(ValueError):
        aggregate_tvar(_g(), 0.0)
