"""Skewness- and kurtosis-risk premia: realized versus BKM-implied moments.

Just as the variance risk premium is implied-minus-realized variance, the
higher-moment risk premia compare the risk-neutral skewness/kurtosis the option
smile implies (Bakshi-Kapadia-Madan, :mod:`quantforge.bkm`) with what actually
realized in the price history (:mod:`quantforge.gramcharlier`). A more negative
risk-neutral skew than realized -- the usual case -- is the market paying up for
crash protection.

This module wires the two together into a single report. Pure standard library.
"""

import math
from typing import Sequence

from .bkm import bkm_moments_from_smile
from .gramcharlier import realized_skewness, realized_excess_kurtosis


def moment_risk_premia(closes: Sequence[float], S0, t, r, vol_fn, q=0.0,
                       n_strikes=401, width=8.0):
    """Realized vs risk-neutral (BKM) skewness and excess kurtosis.

    Args:
        closes: realized price history over the measurement window (its log
            returns give the realized moments).
        S0, t, r, vol_fn, q: inputs for the BKM risk-neutral moments implied by
            the option smile at horizon ``t``.

    Returns a dict with ``realized_skew``, ``implied_skew``, ``skew_premium``
    (implied - realized), and the analogous ``*_kurt`` excess-kurtosis fields.
    """
    closes = [float(c) for c in closes]
    if len(closes) < 4:
        raise ValueError("need at least four closes")
    rets = [math.log(closes[i] / closes[i - 1]) for i in range(1, len(closes))]

    rskew = realized_skewness(rets)
    rkurt = realized_excess_kurtosis(rets)
    _var, iskew, ikurt = bkm_moments_from_smile(S0, t, r, vol_fn, q=q,
                                                n_strikes=n_strikes, width=width)
    return {
        "realized_skew": rskew,
        "implied_skew": iskew,
        "skew_premium": iskew - rskew,
        "realized_kurt": rkurt,
        "implied_kurt": ikurt,
        "kurt_premium": ikurt - rkurt,
    }
