"""Variance risk premium: realized versus implied variance.

The variance risk premium (VRP) is the gap between the variance the market
*implied* (the fair variance-swap strike, a risk-neutral expectation) and the
variance that was actually *realized*. Empirically implied variance exceeds
realized on average -- variance-swap buyers pay a premium for the crash
insurance -- so the realized-minus-implied VRP is typically negative (and the
short-variance carry positive).

This module computes the realized variance from a price history and combines it
with a supplied implied variance (e.g. from :func:`quantforge.variance_swap_from_smile`)
into the additive premium and its ratio. Pure standard library.
"""

import math
from typing import Sequence


def realized_variance(closes: Sequence[float], periods_per_year: int = 252):
    """Annualized realized variance of log returns from a close series."""
    closes = [float(c) for c in closes]
    if len(closes) < 3:
        raise ValueError("need at least three closes")
    rets = [math.log(closes[i] / closes[i - 1]) for i in range(1, len(closes))]
    n = len(rets)
    mean = sum(rets) / n
    var = sum((x - mean) ** 2 for x in rets) / (n - 1)
    return var * periods_per_year


def variance_risk_premium(closes, implied_variance, periods_per_year=252):
    """Realized-minus-implied variance risk premium and its components.

    Args:
        closes: realized price history over the measurement window.
        implied_variance: the annualized implied (variance-swap) variance for
            the same horizon -- e.g. from a smile replication.
        periods_per_year: sampling frequency of the closes.

    Returns a dict with ``realized_variance``, ``implied_variance``,
    ``vrp`` (realized - implied; usually negative), ``ratio``
    (realized / implied), and ``vol_premium`` (implied vol - realized vol, the
    usual positive number quoted in vol points).
    """
    rv = realized_variance(closes, periods_per_year)
    iv = float(implied_variance)
    if iv <= 0:
        raise ValueError("implied variance must be positive")
    return {
        "realized_variance": rv,
        "implied_variance": iv,
        "vrp": rv - iv,
        "ratio": rv / iv,
        "vol_premium": math.sqrt(iv) - math.sqrt(max(rv, 0.0)),
    }
