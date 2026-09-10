"""Ho-Lee (1986) short-rate model: zero-coupon bonds.

Ho-Lee is the simplest no-arbitrage short-rate model: the rate has a
time-dependent drift and constant volatility, with no mean reversion:

    dr = theta(t) dt + sigma dW.

With a flat initial short rate ``r0`` (the constant-drift special case used
here, ``theta(t) = theta``), the zero-coupon bond has the affine closed form

    P(0, t) = exp(-r0 t - 0.5 theta t^2 + sigma^2 t^3 / 6),

so ``B(t) = t`` (the rate has unit sensitivity per year of maturity) and the
volatility contributes a ``+sigma^2 t^3 / 6`` convexity lift. Pure standard
library.
"""

import math


def holee_zero_coupon_bond(r0, t, theta, sigma):
    """Ho-Lee zero-coupon bond price P(0, t) with constant drift ``theta``.

    ``P = exp(-r0 t - 0.5 theta t^2 + sigma^2 t^3 / 6)``.
    """
    if t < 0:
        raise ValueError("t must be non-negative")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if t == 0:
        return 1.0
    return math.exp(-r0 * t - 0.5 * theta * t * t + sigma * sigma * t ** 3 / 6.0)


def holee_zero_coupon_yield(r0, t, theta, sigma):
    """Continuously-compounded yield of the Ho-Lee zero-coupon bond to ``t``.

    ``y(t) = r0 + 0.5 theta t - sigma^2 t^2 / 6`` (linear-in-t drift, quadratic
    convexity pull-down).
    """
    if t <= 0:
        raise ValueError("t must be positive")
    return r0 + 0.5 * theta * t - sigma * sigma * t * t / 6.0
