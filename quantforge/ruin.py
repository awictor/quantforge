"""Risk-of-ruin and drawdown-probability formulas for a trading account.

Two classical results bound the chance an account is wiped out or suffers a given
drawdown before it grows without bound.

The gambler's-ruin problem: starting with ``i`` units and betting one unit per
round with win probability ``p``, the probability of reaching ``0`` before a target
``N`` has the closed form

    q = (1 - p) / p,
    P_ruin = (q^N - q^i) / (q^N - 1)   (q != 1),   or   1 - i/N   (fair game).

For a continuously-compounded account whose log-value follows a Brownian motion
with drift ``mu`` and volatility ``sigma`` (``mu > 0``), the probability that the
equity ever falls to a fraction ``1 - d`` of its starting value -- a loss of at
least the fraction ``d`` at some future time -- is the first-passage law

    P(equity ever <= (1 - d) * start) = (1 - d)^{2 mu / sigma^2}.

This follows from the running-minimum distribution of a Brownian motion with
positive drift, ``P(inf_t X_t <= -a) = exp(-2 mu a / sigma^2)`` with the log-loss
level ``a = -ln(1 - d)``. Pure standard library.
"""

import math


def gamblers_ruin_probability(start, target, win_prob):
    """Probability of hitting 0 before ``target`` in a unit-stake random walk.

    Starting from integer wealth ``start`` with per-round win probability
    ``win_prob`` (win +1, lose -1) and an absorbing target ``target``. The classic
    gambler's-ruin formula; for a fair game (``win_prob = 0.5``) it is the linear
    ``1 - start/target``. Ruin is certain against an unfavorable game as the target
    grows.
    """
    if not (0 < start < target):
        raise ValueError("require 0 < start < target")
    if not (0.0 < win_prob < 1.0):
        raise ValueError("win_prob must be in (0, 1)")
    if win_prob == 0.5:
        return 1.0 - start / target
    q = (1.0 - win_prob) / win_prob
    return (q ** target - q ** start) / (q ** target - 1.0)


def risk_of_ruin_units(win_prob, target_units):
    """Risk of ruin starting with ``target_units`` units against a one-unit target run.

    A trader with a bankroll of ``target_units`` betting one unit at a time with
    edge ``win_prob``: the probability of losing the whole bankroll before doubling
    it is :func:`gamblers_ruin_probability` from ``target_units`` toward
    ``2 * target_units``. A convenient bankroll-in-units risk measure; falls fast in
    the number of units when the game is favorable.
    """
    if not (0.0 < win_prob < 1.0):
        raise ValueError("win_prob must be in (0, 1)")
    if target_units < 1:
        raise ValueError("target_units must be at least 1")
    return gamblers_ruin_probability(target_units, 2 * target_units, win_prob)


def ruin_probability_gbm(loss, mu, sigma):
    """Probability a drifting log-equity ever falls by at least the fraction ``loss``.

    For log-equity following a Brownian motion with drift ``mu > 0`` and volatility
    ``sigma``, the chance the account ever drops to ``(1 - loss)`` of its starting
    value is the first-passage law ``(1 - loss)^{2 mu / sigma^2}``. ``loss`` is a
    fraction in ``(0, 1)``. Decreasing in the drift-to-variance ratio: a stronger
    edge makes a given loss less likely. Returns 1 for a non-positive drift (a
    driftless or losing account eventually hits any loss level almost surely).
    """
    if not (0.0 < loss < 1.0):
        raise ValueError("loss must be in (0, 1)")
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    if mu <= 0.0:
        return 1.0
    exponent = 2.0 * mu / (sigma * sigma)
    return (1.0 - loss) ** exponent
