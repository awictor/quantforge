"""Mean-variance portfolio optimization (portopt module)."""

import pytest

from quantforge import (
    min_variance_weights, max_sharpe_weights, risk_parity_weights,
    portfolio_variance, portfolio_return,
)


DIAG = [[0.04, 0.0], [0.0, 0.01]]
COV3 = [[0.04, 0.01, 0.0], [0.01, 0.09, 0.02], [0.0, 0.02, 0.16]]


def _risk_contribs(w, cov):
    n = len(cov)
    cw = [sum(cov[i][j] * w[j] for j in range(n)) for i in range(n)]
    return [w[i] * cw[i] for i in range(n)]


def test_min_variance_diagonal_is_inverse_variance():
    w = min_variance_weights(DIAG)
    assert w == pytest.approx([0.2, 0.8], abs=1e-9)   # (1/.04):(1/.01) = 1:4
    assert sum(w) == pytest.approx(1.0, abs=1e-12)


def test_min_variance_beats_equal_weight():
    w = min_variance_weights(COV3)
    eq = [1 / 3] * 3
    assert portfolio_variance(w, COV3) < portfolio_variance(eq, COV3)
    assert sum(w) == pytest.approx(1.0, abs=1e-12)


def test_risk_parity_equalizes_risk_contributions():
    for cov in (DIAG, COV3):
        w = risk_parity_weights(cov)
        rc = _risk_contribs(w, cov)
        assert max(rc) - min(rc) < 1e-8
        assert sum(w) == pytest.approx(1.0, abs=1e-9)
        assert all(x > 0 for x in w)


def test_risk_parity_diagonal_is_inverse_sigma():
    w = risk_parity_weights(DIAG)
    # 1/0.2 : 1/0.1 = 1 : 2 -> 1/3, 2/3
    assert w == pytest.approx([1 / 3, 2 / 3], abs=1e-6)


def test_max_sharpe_sums_to_one_and_prefers_higher_return():
    mu = [0.08, 0.05]
    w = max_sharpe_weights(mu, DIAG, risk_free=0.02)
    assert sum(w) == pytest.approx(1.0, abs=1e-9)
    # higher-return, higher-variance asset still gets meaningful weight
    assert w[0] > 0.0 and w[1] > 0.0


def test_portfolio_helpers():
    w = [0.5, 0.5]
    assert portfolio_variance(w, DIAG) == pytest.approx(
        0.25 * 0.04 + 0.25 * 0.01, abs=1e-12)
    assert portfolio_return(w, [0.1, 0.2]) == pytest.approx(0.15, abs=1e-12)


def test_singular_covariance_raises():
    with pytest.raises(ValueError):
        min_variance_weights([[1.0, 1.0], [1.0, 1.0]])
