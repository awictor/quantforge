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


def coupon_bond_option(r0, t_option, cashflows, strike, kappa, theta, sigma,
                       option_type=OptionType.CALL):
    """European option on a coupon bond under Vasicek (Jamshidian decomposition).

    ``cashflows`` is a list of ``(t_i, c_i)`` pairs with ``t_i > t_option``: the
    underlying coupon bond pays ``c_i`` at each ``t_i`` (the last usually
    includes the principal). The option pays ``max(B(t_option) - strike, 0)``
    (call) on the bond's value ``B``.

    Since the Vasicek short rate is one-factor and every zero-coupon bond is
    monotone decreasing in ``r``, Jamshidian's trick applies: find the critical
    rate ``r*`` where the bond value at expiry equals ``strike``, split ``strike``
    into per-cashflow strikes ``K_i = P(t_option, t_i | r*)``, and the coupon-bond
    option is the ``c_i``-weighted sum of zero-coupon-bond options struck at each
    ``K_i``. Exact (no simulation).
    """
    ot = _coerce_type(option_type)
    cfs = sorted(cashflows)
    if not cfs or any(ti <= t_option for ti, _ in cfs):
        raise ValueError("all cashflow times must exceed t_option")

    def bond_value_at(r):
        return sum(c * zero_coupon_bond(r, ti - t_option, kappa, theta, sigma)
                   for ti, c in cfs)

    # Bisection for r* : bond_value_at(r*) = strike (value is decreasing in r).
    lo, hi = -1.0, 1.0
    # Expand until the strike is bracketed (value decreasing).
    while bond_value_at(lo) < strike:
        lo -= 1.0
        if lo < -50.0:
            break
    while bond_value_at(hi) > strike:
        hi += 1.0
        if hi > 50.0:
            break
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if bond_value_at(mid) > strike:
            lo = mid
        else:
            hi = mid
    rstar = 0.5 * (lo + hi)

    total = 0.0
    for ti, c in cfs:
        Ki = zero_coupon_bond(rstar, ti - t_option, kappa, theta, sigma)
        total += c * bond_option(r0, t_option, ti, Ki, kappa, theta, sigma, ot)
    return total


def bond_option_greeks(r0, t_option, t_bond, strike, kappa, theta, sigma,
                       option_type=OptionType.CALL):
    """Greeks of a Vasicek zero-coupon-bond option by central finite differences.

    Central differences of :func:`bond_option` for the short-rate sensitivities
    ``rho_r`` (dV/dr0) and ``gamma_r`` (d2V/dr0^2), and the vol sensitivity
    ``vega`` (dV/dsigma). A bond call *falls* as the short rate rises (higher
    rates discount the bond harder), so ``rho_r < 0`` for a call. Returns a dict
    with ``price``, ``rho_r``, ``gamma_r``, ``vega``.
    """
    ot = _coerce_type(option_type)
    if not (0 < t_option < t_bond):
        raise ValueError("require 0 < t_option < t_bond")

    def px(r=r0, sig=sigma):
        return bond_option(r, t_option, t_bond, strike, kappa, theta, sig, ot)

    base = px()
    hr = 1e-5
    up, dn = px(r=r0 + hr), px(r=r0 - hr)
    rho_r = (up - dn) / (2.0 * hr)
    gamma_r = (up - 2.0 * base + dn) / (hr * hr)
    hv = 1e-6
    vega = (px(sig=sigma + hv) - px(sig=max(sigma - hv, 0.0))) / (
        (2.0 * hv) if sigma - hv >= 0 else hv)
    return {"price": base, "rho_r": rho_r, "gamma_r": gamma_r, "vega": vega}


def bond_greeks(r0, t, kappa, theta, sigma):
    """Exact rate sensitivities of a Vasicek zero-coupon bond.

    ``P = A(t) e^{-B(t) r0}`` (with ``B = _B(kappa, t)``), so ``rho_r = -B P``,
    ``gamma_r = B^2 P``, rate ``duration = B``, and ``convexity = B^2``. Returns
    a dict with ``price``, ``rho_r``, ``gamma_r``, ``duration``, ``convexity``.
    """
    if t < 0:
        raise ValueError("t must be non-negative")
    price = zero_coupon_bond(r0, t, kappa, theta, sigma)
    if t == 0:
        return {"price": 1.0, "rho_r": 0.0, "gamma_r": 0.0,
                "duration": 0.0, "convexity": 0.0}
    B = _B(kappa, t)
    return {"price": price, "rho_r": -B * price, "gamma_r": B * B * price,
            "duration": B, "convexity": B * B}
