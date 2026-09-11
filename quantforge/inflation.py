"""Inflation-linked bond and breakeven-inflation analytics.

Real (TIPS-style) instruments scale their cashflows by the ratio of the
reference price index at settlement to the index at issue. This module handles
the index ratio, the inflation-adjusted principal, the Fisher relation between
nominal and real rates, and the breakeven inflation implied by a nominal/real
yield pair. Pure standard library.
"""

import math


def index_ratio(index_settle, index_base) -> float:
    """Index ratio ``CPI_settle / CPI_base`` used to inflate the principal.

    Above 1 when the price index has risen since issue. Both indices must be
    positive.
    """
    if index_base <= 0 or index_settle <= 0:
        raise ValueError("index levels must be positive")
    return index_settle / index_base


def inflation_adjusted_principal(face, index_settle, index_base) -> float:
    """Inflation-adjusted principal ``face * index_ratio`` (the linker notional)."""
    return face * index_ratio(index_settle, index_base)


def fisher_real_rate(nominal, inflation) -> float:
    """Exact Fisher real rate ``(1 + nominal)/(1 + inflation) - 1``.

    The rate that, compounded with inflation, reproduces the nominal rate. For
    small rates it is approximately ``nominal - inflation``.
    """
    return (1.0 + nominal) / (1.0 + inflation) - 1.0


def fisher_nominal_rate(real, inflation) -> float:
    """Exact Fisher nominal rate ``(1 + real)(1 + inflation) - 1``."""
    return (1.0 + real) * (1.0 + inflation) - 1.0


def breakeven_inflation(nominal_yield, real_yield) -> float:
    """Breakeven inflation implied by a nominal and a real yield (Fisher).

    ``(1 + nominal)/(1 + real) - 1`` -- the inflation rate at which a nominal and
    an inflation-linked bond of the same maturity have equal return. The market's
    inflation expectation (plus risk premium).
    """
    if real_yield <= -1.0:
        raise ValueError("real_yield must exceed -100%")
    return (1.0 + nominal_yield) / (1.0 + real_yield) - 1.0


def real_from_breakeven(nominal_yield, breakeven) -> float:
    """Real yield implied by a nominal yield and a breakeven inflation rate.

    ``(1 + nominal)/(1 + breakeven) - 1`` -- inverse of
    :func:`breakeven_inflation`.
    """
    if breakeven <= -1.0:
        raise ValueError("breakeven must exceed -100%")
    return (1.0 + nominal_yield) / (1.0 + breakeven) - 1.0


def linker_price(real_cashflows, real_yield, index_settle, index_base) -> float:
    """Dirty price of an inflation-linked bond off real cashflows.

    ``real_cashflows`` is ``[(t, real_amount), ...]`` in constant (issue-date)
    money. Each flow is discounted at the continuously-compounded ``real_yield``
    and then the whole bond is inflated by the settlement index ratio:

        price = (index_settle / index_base) * sum_i real_amount_i e^{-r t_i}

    Because the index ratio multiplies every flow, the price is degree-one
    homogeneous in it -- stripping the ratio recovers a standard real-yield bond
    price (the invariant tested against :mod:`quantforge.bondmath`).
    """
    ratio = index_ratio(index_settle, index_base)
    pv = sum(amt * math.exp(-real_yield * t) for t, amt in real_cashflows)
    return ratio * pv


def linker_real_yield(real_cashflows, price, index_settle, index_base,
                      tol=1e-10, max_iter=100) -> float:
    """Continuously-compounded real yield reproducing a linker ``price``.

    Deflates the quoted price by the index ratio and solves the standard real-
    cashflow bond yield by bisection (price is monotone decreasing in the yield).
    Inverse of :func:`linker_price`.
    """
    if price <= 0:
        raise ValueError("price must be positive")
    ratio = index_ratio(index_settle, index_base)
    target = price / ratio  # real (deflated) price

    def real_px(r):
        return sum(amt * math.exp(-r * t) for t, amt in real_cashflows)

    lo, hi = -0.5, 5.0
    p_lo, p_hi = real_px(lo), real_px(hi)
    if not (min(p_lo, p_hi) - 1e-9 <= target <= max(p_lo, p_hi) + 1e-9):
        raise ValueError("price outside the achievable yield range")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        pm = real_px(mid)
        if abs(pm - target) < tol:
            return mid
        if pm > target:   # price too high -> raise yield
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
