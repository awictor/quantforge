"""Black-Litterman posterior returns (portopt module)."""

import pytest

from quantforge import (
    implied_equilibrium_returns, black_litterman_returns,
)


COV = [[0.04, 0.01, 0.0], [0.01, 0.09, 0.02], [0.0, 0.02, 0.16]]
MW = [0.5, 0.3, 0.2]


def test_equilibrium_is_lambda_cw():
    pi = implied_equilibrium_returns(COV, MW, risk_aversion=2.5)
    cw = [sum(COV[i][j] * MW[j] for j in range(3)) for i in range(3)]
    assert pi == pytest.approx([2.5 * v for v in cw], abs=1e-12)


def test_no_views_returns_prior():
    pi = implied_equilibrium_returns(COV, MW)
    post = black_litterman_returns(COV, MW, [], [])
    assert post == pytest.approx(pi, abs=1e-12)


def test_bullish_view_raises_that_asset():
    pi = implied_equilibrium_returns(COV, MW)
    P = [[1.0, 0.0, 0.0]]
    Q = [pi[0] * 1.5]   # view above the prior
    post = black_litterman_returns(COV, MW, P, Q)
    assert post[0] > pi[0]


def test_relative_view_moves_spread():
    pi = implied_equilibrium_returns(COV, MW)
    prior_spread = pi[0] - pi[1]
    P = [[1.0, -1.0, 0.0]]
    Q = [0.02]   # asset 0 outperforms asset 1 by 2%
    post = black_litterman_returns(COV, MW, P, Q)
    post_spread = post[0] - post[1]
    # posterior spread lies between the prior spread and the view.
    assert min(prior_spread, 0.02) <= post_spread <= max(prior_spread, 0.02)


def test_higher_confidence_view_pulls_closer():
    # A tighter omega (smaller variance) pulls the posterior nearer the view.
    P = [[1.0, 0.0, 0.0]]
    Q = [0.15]
    loose = black_litterman_returns(COV, MW, P, Q, omega=[[0.1]])
    tight = black_litterman_returns(COV, MW, P, Q, omega=[[1e-6]])
    assert abs(tight[0] - 0.15) < abs(loose[0] - 0.15)


def test_validation():
    with pytest.raises(ValueError):
        black_litterman_returns(COV, MW, [[1.0, 0.0, 0.0]], [0.1, 0.2])  # Q len
    with pytest.raises(ValueError):
        implied_equilibrium_returns(COV, [0.5, 0.5])  # weight length
