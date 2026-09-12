"""Futures/forward convexity adjustment for interest-rate futures.

An interest-rate futures contract is marked to market daily, so its rate exceeds
the forward rate for the same period: a short-rate rise lowers the futures price
exactly when financing the variation margin is cheap, a correlation that makes the
futures rate biased high. Converting a futures rate to a forward rate subtracts a
convexity adjustment. Under Ho-Lee (constant volatility ``sigma``) the adjustment
to a rate fixing at ``t1`` and applying to ``t2`` is

    adj = 0.5 * sigma^2 * t1 * t2,     forward = futures - adj

and under Hull-White (mean reversion ``a``) it uses the exponential terms. Pure
standard library.
"""

import math


def ho_lee_convexity_adjustment(sigma, t1, t2):
    """Ho-Lee convexity adjustment ``0.5 sigma^2 t1 t2`` (continuous-comp rates).

    ``t1`` is the time to the rate fixing, ``t2`` the time to the end of the
    underlying accrual period (``t2 >= t1``). Non-negative, zero at zero vol, and
    growing with both horizons and the volatility.
    """
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if t1 < 0 or t2 < t1:
        raise ValueError("require 0 <= t1 <= t2")
    return 0.5 * sigma * sigma * t1 * t2


def hull_white_convexity_adjustment(sigma, a, t1, t2):
    """Hull-White (mean-reverting) convexity adjustment.

    With mean-reversion speed ``a`` the adjustment is

        adj = (sigma^2 / (2 a^2)) * (1 - e^{-a(t2 - t1)}) *
              [ (1 - e^{-a(t2 - t1)}) * (1 - e^{-2 a t1}) / (2 a)
                + a * B(0, t1)^2 ... ]  (standard HW futures-forward formula)

    Implemented in the common compact form
        adj = (sigma^2 / (2 a)) * B(t1, t2) * [ B(t1, t2) (1 - e^{-2 a t1})
              + 2 a B(0, t1)^2 ] / 2,
    with ``B(u, v) = (1 - e^{-a (v - u)}) / a``. Reduces to
    :func:`ho_lee_convexity_adjustment` as ``a -> 0``.
    """
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if t1 < 0 or t2 < t1:
        raise ValueError("require 0 <= t1 <= t2")
    if a < 0:
        raise ValueError("a must be non-negative")
    if abs(a) < 1e-8:
        return ho_lee_convexity_adjustment(sigma, t1, t2)
    B12 = (1.0 - math.exp(-a * (t2 - t1))) / a
    B01 = (1.0 - math.exp(-a * t1)) / a
    term = B12 * (1.0 - math.exp(-2.0 * a * t1)) + 2.0 * a * B01 * B01
    return (sigma * sigma / (2.0 * a)) * B12 * term / 2.0


def forward_from_futures(futures_rate, sigma, t1, t2, a=0.0):
    """Forward rate from a futures rate, subtracting the convexity adjustment.

    Uses Hull-White (or Ho-Lee when ``a = 0``). The forward is below the futures
    rate by the (non-negative) adjustment.
    """
    if a and a > 1e-8:
        adj = hull_white_convexity_adjustment(sigma, a, t1, t2)
    else:
        adj = ho_lee_convexity_adjustment(sigma, t1, t2)
    return futures_rate - adj


def futures_from_forward(forward_rate, sigma, t1, t2, a=0.0):
    """Futures rate from a forward rate, adding the convexity adjustment.

    Inverse of :func:`forward_from_futures`; the futures rate is above the forward.
    """
    if a and a > 1e-8:
        adj = hull_white_convexity_adjustment(sigma, a, t1, t2)
    else:
        adj = ho_lee_convexity_adjustment(sigma, t1, t2)
    return forward_rate + adj
