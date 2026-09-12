"""Generalized downside-risk ratios: Kappa and the upside-potential ratio.

Downside-risk performance measures judge return per unit of *shortfall* below a
target ``tau``, not per unit of total volatility. The lower partial moment of
order ``n`` is

    LPM_n(tau) = mean( max(tau - r, 0)^n ).

Kaplan-Knowles' Kappa of order ``n`` normalizes the excess return by its ``n``-th
root:

    kappa_n = (mean(r) - tau) / LPM_n(tau)^(1/n).

``kappa_1`` is the Omega-Sharpe ratio and ``kappa_2`` is the (per-period) Sortino
ratio, so Kappa unifies them; higher ``n`` weights deep shortfalls more heavily.
The upside-potential ratio replaces the numerator with the *upside* expectation,
``mean(max(r - tau, 0))``, over the downside deviation -- reward for the chance of
beating the target per unit of downside. Pure standard library.
"""

import math


def _mean(x):
    return sum(x) / len(x)


def lower_partial_moment(returns, tau=0.0, order=2):
    """``n``-th lower partial moment about ``tau``: mean(max(tau - r, 0)^n)."""
    if len(returns) == 0:
        raise ValueError("need at least one return")
    if order <= 0:
        raise ValueError("order must be positive")
    return sum(max(tau - r, 0.0) ** order for r in returns) / len(returns)


def kappa_ratio(returns, tau=0.0, order=2):
    """Kaplan-Knowles Kappa of the given ``order`` about target ``tau``.

    ``kappa_n = (mean - tau) / LPM_n^(1/n)``. ``order=1`` is the Omega-Sharpe
    ratio, ``order=2`` the per-period Sortino ratio. Raises if there is no
    downside (zero lower partial moment).
    """
    if len(returns) < 1:
        raise ValueError("need at least one return")
    lpm = lower_partial_moment(returns, tau, order)
    if lpm <= 0.0:
        raise ValueError("no downside below tau; kappa undefined")
    return (_mean(returns) - tau) / (lpm ** (1.0 / order))


def upside_potential_ratio(returns, tau=0.0):
    """Upside-potential ratio: mean upside over downside deviation about ``tau``.

    ``UPR = mean(max(r - tau, 0)) / sqrt(LPM_2(tau))`` -- expected outperformance
    of the target per unit of downside risk. Raises if there is no downside.
    """
    if len(returns) < 1:
        raise ValueError("need at least one return")
    upside = sum(max(r - tau, 0.0) for r in returns) / len(returns)
    dd = math.sqrt(lower_partial_moment(returns, tau, 2))
    if dd <= 0.0:
        raise ValueError("no downside below tau; ratio undefined")
    return upside / dd
