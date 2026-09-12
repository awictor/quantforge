"""FX forwards under covered interest parity (CIP).

Quotes are ``price`` currency per unit of ``base`` currency (e.g. USD per EUR
for EURUSD). With continuously-compounded domestic (price-currency) rate
``r_price`` and foreign (base-currency) rate ``r_base`` over ``t`` years, no-
arbitrage forces the forward

    F = S * exp((r_price - r_base) * t).

Forward points are ``F - S``; an FX swap exchanges spot for forward at those
points. Inverting CIP recovers an implied rate from a quoted forward. Pure
standard library.
"""

import math


def cross_rate(ab, cb, via_is_quote=True):
    """Cross rate from two rates sharing a common currency.

    Given ``A/B`` and ``C/B`` (both quoted against the common currency ``B``), the
    cross ``A/C = (A/B) / (C/B)`` when ``via_is_quote`` (the common currency is the
    quote of both), else ``A/C = (A/B) * (B/C)`` for ``A/B`` and ``B/C``. Positive
    rates required.
    """
    if ab <= 0 or cb <= 0:
        raise ValueError("rates must be positive")
    return ab / cb if via_is_quote else ab * cb


def triangular_arbitrage(ab, bc, ca):
    """Triangular-arbitrage profit factor around a currency loop ``A->B->C->A``.

    Converting one unit of A through ``A/B``... actually multiplying the three
    quoted legs ``(A per B) (B per C) (C per A)`` returns the units of A after a
    round trip; it equals one in an arbitrage-free market. Returns the product;
    values above one (net of costs) signal a profitable loop, below one the reverse
    direction.
    """
    if ab <= 0 or bc <= 0 or ca <= 0:
        raise ValueError("rates must be positive")
    return ab * bc * ca


def is_arbitrage_free(ab, bc, ca, tol=1e-9):
    """True if the triangular loop is arbitrage-free within ``tol``.

    Checks ``|triangular_arbitrage - 1| <= tol``; the cross rates are mutually
    consistent when the round-trip product is one.
    """
    return abs(triangular_arbitrage(ab, bc, ca) - 1.0) <= tol


def fx_forward(spot, r_price, r_base, t) -> float:
    """Covered-interest-parity forward FX rate ``S exp((r_price - r_base) t)``.

    A base currency yielding more than the price currency (``r_base > r_price``)
    trades at a forward discount (``F < S``), and vice versa.
    """
    if spot <= 0:
        raise ValueError("spot must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    return spot * math.exp((r_price - r_base) * t)


def fx_forward_from_curves(spot, price_curve, base_curve, t) -> float:
    """FX forward from two discount curves: ``S * DF_base(t) / DF_price(t)``.

    Curve-based covered interest parity: a base currency that discounts more
    steeply (higher rates, lower ``DF_base``) trades at a forward discount. Each
    curve is anything callable as ``curve.df(t)`` (e.g.
    :class:`quantforge.DiscountCurve`) or a plain ``curve(t)`` returning the
    discount factor. Consistent with :func:`fx_forward` when the curves are flat
    exponentials.
    """
    if spot <= 0:
        raise ValueError("spot must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    dfp = price_curve.df if hasattr(price_curve, "df") else price_curve
    dfb = base_curve.df if hasattr(base_curve, "df") else base_curve
    dp = dfp(t)
    if dp <= 0.0:
        raise ValueError("price-curve discount factor must be positive")
    return spot * dfb(t) / dp


def forward_points(spot, r_price, r_base, t) -> float:
    """Forward points ``F - S`` (positive when the base is at a forward premium)."""
    return fx_forward(spot, r_price, r_base, t) - spot


def fx_swap_points(spot, r_price, r_base, t_near, t_far) -> float:
    """FX-swap points between two tenors: ``F(t_far) - F(t_near)``.

    The pips exchanged in a forward-forward FX swap rolling from the near to the
    far date.
    """
    if not (0.0 <= t_near < t_far):
        raise ValueError("require 0 <= t_near < t_far")
    return (fx_forward(spot, r_price, r_base, t_far)
            - fx_forward(spot, r_price, r_base, t_near))


def implied_base_rate(spot, forward, r_price, t) -> float:
    """Base-currency rate implied by a quoted forward (invert CIP).

    ``r_base = r_price - ln(forward / spot) / t``.
    """
    if spot <= 0 or forward <= 0:
        raise ValueError("spot and forward must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    return r_price - math.log(forward / spot) / t


def implied_price_rate(spot, forward, r_base, t) -> float:
    """Price-currency rate implied by a quoted forward (invert CIP).

    ``r_price = r_base + ln(forward / spot) / t``.
    """
    if spot <= 0 or forward <= 0:
        raise ValueError("spot and forward must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    return r_base + math.log(forward / spot) / t
