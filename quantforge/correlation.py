"""Implied correlation of an index from its members' and the index's vols.

An index's variance decomposes into its constituents' variances plus every
pairwise covariance:

    sigma_I^2 = sum_i w_i^2 sigma_i^2
              + sum_{i != j} w_i w_j sigma_i sigma_j rho_{ij}.

Assuming a single common correlation ``rho`` for every pair, this inverts to a
closed form for the **implied correlation** consistent with a quoted index vol:

    rho = (sigma_I^2 - sum_i w_i^2 sigma_i^2)
          / (sum_{i != j} w_i w_j sigma_i sigma_j).

The denominator is ``(sum_i w_i sigma_i)^2 - sum_i w_i^2 sigma_i^2``. This is
the standard dispersion-trading measure: high implied correlation means the
index is expensive relative to its members (or vice-versa).
"""

import math
from typing import Sequence


def _validate(weights, vols):
    n = len(weights)
    if n != len(vols) or n < 2:
        raise ValueError("need at least two matching (weight, vol) entries")
    for v in vols:
        if v < 0:
            raise ValueError("vols must be non-negative")


def index_vol_from_correlation(weights: Sequence[float], vols: Sequence[float],
                               rho: float) -> float:
    """Index volatility implied by member weights/vols and a common correlation."""
    _validate(weights, vols)
    if not (-1.0 <= rho <= 1.0):
        raise ValueError("rho must be in [-1, 1]")
    diag = sum(w * w * s * s for w, s in zip(weights, vols))
    wsum = sum(w * s for w, s in zip(weights, vols))
    cross = wsum * wsum - diag           # sum_{i!=j} w_i w_j s_i s_j
    var = diag + rho * cross
    return math.sqrt(max(var, 0.0))


def implied_correlation(weights: Sequence[float], vols: Sequence[float],
                        index_vol: float) -> float:
    """Common implied correlation consistent with the quoted ``index_vol``.

    Returns rho in principle within [-1, 1]; a value outside that band signals
    an index vol inconsistent with the member vols (arbitrage or stale quotes)
    and is returned unclamped so the caller can see it.
    """
    _validate(weights, vols)
    if index_vol < 0:
        raise ValueError("index_vol must be non-negative")
    diag = sum(w * w * s * s for w, s in zip(weights, vols))
    wsum = sum(w * s for w, s in zip(weights, vols))
    cross = wsum * wsum - diag
    if abs(cross) < 1e-15:
        raise ValueError("degenerate basket (zero cross term); cannot imply rho")
    return (index_vol * index_vol - diag) / cross


def dispersion_basket_vol(weights: Sequence[float], vols: Sequence[float]) -> float:
    """The zero-correlation ("fully diversified") index vol, sqrt(sum w^2 sig^2).

    A useful lower reference: the index vol if the members were uncorrelated.
    """
    _validate(weights, vols)
    return math.sqrt(sum(w * w * s * s for w, s in zip(weights, vols)))
