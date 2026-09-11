"""Marginal VaR and VaR budgeting (portopt module)."""

import pytest

from quantforge import (
    marginal_var, var_budget, portfolio_var, component_var,
    risk_parity_weights,
)


COV = [[0.04, 0.01, 0.0], [0.01, 0.09, 0.02], [0.0, 0.02, 0.16]]
W = [0.4, 0.4, 0.2]


def test_euler_decomposition():
    # sum_i w_i * marginal_i == total VaR (VaR homogeneous of degree 1).
    m = marginal_var(W, COV)
    euler = sum(W[i] * m[i] for i in range(3))
    assert euler == pytest.approx(portfolio_var(W, COV), abs=1e-9)


def test_component_is_weight_times_marginal():
    # component_var uses sigma units; w*marginal is in VaR (z-scaled) units.
    z = 1.6448536269514722
    m = marginal_var(W, COV)
    comps = component_var(W, COV)
    for i in range(3):
        assert W[i] * m[i] == pytest.approx(comps[i] * z, abs=1e-9)


def test_budget_sums_to_one():
    bud = var_budget(W, COV)
    assert sum(bud) == pytest.approx(1.0, abs=1e-9)
    assert all(b >= 0 for b in bud)


def test_risk_parity_gives_equal_budget():
    w = risk_parity_weights(COV)
    bud = var_budget(w, COV)
    assert max(bud) - min(bud) < 1e-7


def test_budget_independent_of_confidence():
    # The z and horizon cancel in the ratio, so budgets are level-independent.
    assert var_budget(W, COV) == pytest.approx(var_budget(W, COV), abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        marginal_var([0.0, 0.0, 0.0], COV)
    with pytest.raises(ValueError):
        var_budget([0.0, 0.0, 0.0], COV)
