"""Risk budgeting: contributions and target-budget weights."""

import pytest

from quantforge import (risk_contributions, risk_budget_weights,
                        risk_parity_weights)
from quantforge.portopt import portfolio_variance

COV = [[0.04, 0.006, 0.0], [0.006, 0.02, 0.004], [0.0, 0.004, 0.09]]


def test_equal_budget_is_risk_parity():
    w_eq = risk_budget_weights(COV, [1, 1, 1])
    w_rp = risk_parity_weights(COV)
    assert all(abs(w_eq[i] - w_rp[i]) < 1e-6 for i in range(3))


def test_contributions_match_target():
    budgets = [0.5, 0.3, 0.2]
    w = risk_budget_weights(COV, budgets)
    rc = risk_contributions(w, COV)
    total = sum(rc)
    pct = [r / total for r in rc]
    for i in range(3):
        assert abs(pct[i] - budgets[i]) < 1e-4


def test_contributions_sum_to_variance():
    w = risk_budget_weights(COV, [0.4, 0.4, 0.2])
    assert abs(sum(risk_contributions(w, COV)) - portfolio_variance(w, COV)) < 1e-12


def test_weights_positive_and_sum_one():
    w = risk_budget_weights(COV, [0.6, 0.2, 0.2])
    assert abs(sum(w) - 1.0) < 1e-9
    assert all(x > 0 for x in w)


def test_higher_budget_more_weight_for_equal_vol():
    # Equal variances, no correlation: a larger budget gets proportionally more weight.
    cov = [[0.04, 0.0], [0.0, 0.04]]
    w = risk_budget_weights(cov, [0.75, 0.25])
    assert w[0] > w[1]


def test_validation():
    with pytest.raises(ValueError):
        risk_budget_weights(COV, [1, 1])            # wrong length
    with pytest.raises(ValueError):
        risk_budget_weights(COV, [0.5, 0.5, 0.0])   # non-positive budget
    with pytest.raises(ValueError):
        risk_contributions([0.5, 0.5], COV)         # wrong weight length
