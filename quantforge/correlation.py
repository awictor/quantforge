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


def dispersion_trade_pnl(weights, realized_member_vols, realized_index_vol,
                         strike_member_vols, strike_index_vol,
                         variance_notional=1.0):
    """P&L of a variance dispersion trade (short index var, long member var).

    A dispersion trade sells index variance and buys the weighted member
    variances. Its variance P&L per unit notional is

        (sum_i w_i (rv_i^2 - k_i^2))  -  (rv_index^2 - k_index^2),

    the long member-variance legs minus the short index-variance leg (strikes
    ``k``). Because index variance carries the correlation, the trade profits when
    realized correlation comes in *below* what was implied (index realizes calmer
    than the members would imply), and is zero when realized matches strikes.
    """
    _validate(weights, realized_member_vols)
    if len(strike_member_vols) != len(weights):
        raise ValueError("strike_member_vols must align with weights")
    member_leg = sum(w * (rv * rv - k * k)
                     for w, rv, k in zip(weights, realized_member_vols,
                                         strike_member_vols))
    index_leg = realized_index_vol ** 2 - strike_index_vol ** 2
    return variance_notional * (member_leg - index_leg)


def dispersion_basket_vol(weights: Sequence[float], vols: Sequence[float]) -> float:
    """The zero-correlation ("fully diversified") index vol, sqrt(sum w^2 sig^2).

    A useful lower reference: the index vol if the members were uncorrelated.
    """
    _validate(weights, vols)
    return math.sqrt(sum(w * w * s * s for w, s in zip(weights, vols)))


def correlation_term_structure(weights, member_vol_curves, index_vol_curve,
                               expiries):
    """Implied correlation at each expiry across a term structure.

    Args:
        weights: index member weights (constant across expiries).
        member_vol_curves: list per member of that member's vol at each expiry,
            i.e. member_vol_curves[i][j] is member i's vol at expiries[j].
        index_vol_curve: the index's implied vol at each expiry.
        expiries: the tenors (used only as labels in the returned pairs).

    Returns a list of ``(expiry, implied_correlation)`` pairs, applying
    :func:`implied_correlation` slice by slice.
    """
    n_exp = len(expiries)
    if len(index_vol_curve) != n_exp:
        raise ValueError("index_vol_curve must match expiries length")
    for curve in member_vol_curves:
        if len(curve) != n_exp:
            raise ValueError("each member vol curve must match expiries length")

    out = []
    for j in range(n_exp):
        vols_j = [member_vol_curves[i][j] for i in range(len(weights))]
        rho = implied_correlation(weights, vols_j, index_vol_curve[j])
        out.append((expiries[j], rho))
    return out


def ewma_covariance(returns_x, returns_y, lam=0.94):
    """Exponentially-weighted covariance of two aligned return series.

    RiskMetrics-style recursion ``s_t = lam s_{t-1} + (1-lam) x_t y_t`` seeded
    from the first product, giving more weight to recent observations. ``lam``
    is the decay (0.94 for daily data). The series must be equal length and
    zero-mean is assumed (the RiskMetrics convention for returns).
    """
    if len(returns_x) != len(returns_y):
        raise ValueError("return series must be equal length")
    if len(returns_x) < 2:
        raise ValueError("need at least two observations")
    if not (0.0 < lam < 1.0):
        raise ValueError("lam must be in (0, 1)")
    s = returns_x[0] * returns_y[0]
    for x, y in zip(returns_x[1:], returns_y[1:]):
        s = lam * s + (1.0 - lam) * x * y
    return s


def ewma_correlation(returns_x, returns_y, lam=0.94):
    """Exponentially-weighted correlation of two aligned return series.

    The EWMA covariance divided by the product of the EWMA volatilities (all on
    the same decay), so it stays in ``[-1, 1]``.
    """
    cov = ewma_covariance(returns_x, returns_y, lam)
    vx = ewma_covariance(returns_x, returns_x, lam)
    vy = ewma_covariance(returns_y, returns_y, lam)
    denom = math.sqrt(vx * vy)
    if denom <= 0.0:
        raise ValueError("degenerate (zero-variance) series")
    return cov / denom


def realized_beta(asset_returns, market_returns):
    """Realized beta of an asset to the market: ``Cov(a, m) / Var(m)``.

    Ordinary (equal-weight) sample covariance over variance, the slope of a
    regression of asset returns on market returns. Series must be equal length.
    """
    n = len(asset_returns)
    if n != len(market_returns):
        raise ValueError("return series must be equal length")
    if n < 2:
        raise ValueError("need at least two observations")
    ma = sum(asset_returns) / n
    mm = sum(market_returns) / n
    cov = sum((a - ma) * (m - mm) for a, m in zip(asset_returns, market_returns)) / n
    var = sum((m - mm) ** 2 for m in market_returns) / n
    if var <= 0.0:
        raise ValueError("market variance must be positive")
    return cov / var
