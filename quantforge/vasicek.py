"""Vasicek (1977) short-rate model: zero-coupon bonds and bond options.

The short rate follows a mean-reverting Ornstein-Uhlenbeck process:

    dr = kappa (theta - r) dt + sigma dW

with mean-reversion speed ``kappa``, long-run level ``theta``, and volatility
``sigma``. Zero-coupon bonds have the affine closed form ``P(t,T) = A e^{-B r}``,
and European options on a zero-coupon bond have a Jamshidian closed form (a
Black-Scholes-like formula with the bond's forward-price volatility). Everything
is pure standard library.
"""

import math

from .mathfns import norm_cdf
from .bsm import OptionType, _coerce_type


def _B(kappa, tau):
    if abs(kappa) < 1e-12:
        return tau
    return (1.0 - math.exp(-kappa * tau)) / kappa


def zero_coupon_bond(r0, t, kappa, theta, sigma):
    """Vasicek zero-coupon bond price P(0, t) for a unit face, given r(0)=r0.

    ``P = A(t) * exp(-B(t) * r0)`` with the standard affine coefficients.
    """
    if t < 0:
        raise ValueError("t must be non-negative")
    if t == 0:
        return 1.0
    B = _B(kappa, t)
    if abs(kappa) < 1e-12:
        # kappa -> 0 limit of the A-term.
        lnA = (sigma * sigma * t ** 3) / 6.0
    else:
        lnA = ((theta - sigma * sigma / (2.0 * kappa * kappa)) * (B - t)
               - (sigma * sigma) / (4.0 * kappa) * B * B)
    return math.exp(lnA - B * r0)


def zero_coupon_yield(r0, t, kappa, theta, sigma):
    """Continuously-compounded yield of the Vasicek zero-coupon bond to ``t``."""
    if t <= 0:
        raise ValueError("t must be positive")
    return -math.log(zero_coupon_bond(r0, t, kappa, theta, sigma)) / t


def bond_option(r0, t_option, t_bond, strike, kappa, theta, sigma,
                option_type=OptionType.CALL):
    """European option (Jamshidian) on a zero-coupon bond under Vasicek.

    Args:
        t_option: option expiry. t_bond: the underlying bond's maturity
            (``t_bond > t_option``). strike: strike on the bond price.

    A call pays ``max(P(t_option, t_bond) - strike, 0)`` at the option expiry.
    Uses the closed-form bond-price volatility.
    """
    ot = _coerce_type(option_type)
    if not (0 < t_option < t_bond):
        raise ValueError("require 0 < t_option < t_bond")

    P_bond = zero_coupon_bond(r0, t_bond, kappa, theta, sigma)   # to bond maturity
    P_opt = zero_coupon_bond(r0, t_option, kappa, theta, sigma)  # to option expiry

    B = _B(kappa, t_bond - t_option)
    if abs(kappa) < 1e-12:
        sig_p = sigma * B * math.sqrt(t_option)
    else:
        sig_p = sigma * B * math.sqrt((1.0 - math.exp(-2.0 * kappa * t_option))
                                      / (2.0 * kappa))
    if sig_p < 1e-14:
        # No vol: intrinsic on the forward bond price.
        fwd = P_bond / P_opt
        payoff = max(fwd - strike, 0.0) if ot is OptionType.CALL else max(strike - fwd, 0.0)
        return P_opt * payoff

    d1 = (math.log(P_bond / (strike * P_opt)) + 0.5 * sig_p * sig_p) / sig_p
    d2 = d1 - sig_p
    if ot is OptionType.CALL:
        return P_bond * norm_cdf(d1) - strike * P_opt * norm_cdf(d2)
    return strike * P_opt * norm_cdf(-d2) - P_bond * norm_cdf(-d1)
