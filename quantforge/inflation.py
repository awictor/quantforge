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
