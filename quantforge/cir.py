"""Cox-Ingersoll-Ross (1985) short-rate model: zero-coupon bonds and yields.

CIR is an affine short-rate model like Vasicek but with a square-root diffusion
that keeps the rate non-negative:

    dr = kappa (theta - r) dt + sigma sqrt(r) dW

Zero-coupon bonds have the affine closed form ``P(t) = A(t) e^{-B(t) r0}`` with

    gamma = sqrt(kappa^2 + 2 sigma^2)
    B(t)  = 2 (e^{gamma t} - 1) / ((gamma + kappa)(e^{gamma t} - 1) + 2 gamma)
    A(t)  = [ 2 gamma e^{(kappa+gamma) t / 2}
              / ((gamma + kappa)(e^{gamma t} - 1) + 2 gamma) ]^{2 kappa theta / sigma^2}

Pure standard library.
"""

import math

from .cev import noncentral_chisq_cdf
from .bsm import OptionType, _coerce_type


def _AB(t, kappa, theta, sigma):
    g = math.sqrt(kappa * kappa + 2.0 * sigma * sigma)
    eg = math.exp(g * t)
    denom = (g + kappa) * (eg - 1.0) + 2.0 * g
    B = 2.0 * (eg - 1.0) / denom
    A = (2.0 * g * math.exp((kappa + g) * t / 2.0) / denom) ** (
        2.0 * kappa * theta / (sigma * sigma))
    return A, B


def cir_zero_coupon_bond(r0, t, kappa, theta, sigma):
    """CIR zero-coupon bond price P(0, t) for a unit face, given r(0)=r0.

    Requires ``r0 >= 0`` and positive parameters.
    """
    if t < 0:
        raise ValueError("t must be non-negative")
    if r0 < 0:
        raise ValueError("r0 must be non-negative in CIR")
    if kappa <= 0 or theta < 0 or sigma <= 0:
        raise ValueError("require kappa > 0, theta >= 0, sigma > 0")
    if t == 0:
        return 1.0
    A, B = _AB(t, kappa, theta, sigma)
    return A * math.exp(-B * r0)


def cir_zero_coupon_yield(r0, t, kappa, theta, sigma):
    """Continuously-compounded yield of the CIR zero-coupon bond to ``t``."""
    if t <= 0:
        raise ValueError("t must be positive")
    return -math.log(cir_zero_coupon_bond(r0, t, kappa, theta, sigma)) / t


def cir_bond_greeks(r0, t, kappa, theta, sigma):
    """Rate sensitivities of a CIR zero-coupon bond, exact.

    The bond is ``P = A(t) e^{-B(t) r0}``, so its short-rate sensitivities are
    closed form: ``rho_r = dP/dr0 = -B P`` and ``gamma_r = d2P/dr0^2 = B^2 P``.
    The rate ``duration`` is ``-1/P dP/dr0 = B`` and ``convexity`` is
    ``1/P d2P/dr0^2 = B^2``. Returns a dict with ``price``, ``rho_r``,
    ``gamma_r``, ``duration``, ``convexity``.
    """
    if t < 0:
        raise ValueError("t must be non-negative")
    if r0 < 0:
        raise ValueError("r0 must be non-negative in CIR")
    if kappa <= 0 or theta < 0 or sigma <= 0:
        raise ValueError("require kappa > 0, theta >= 0, sigma > 0")
    price = cir_zero_coupon_bond(r0, t, kappa, theta, sigma)
    if t == 0:
        return {"price": 1.0, "rho_r": 0.0, "gamma_r": 0.0,
                "duration": 0.0, "convexity": 0.0}
    _A, B = _AB(t, kappa, theta, sigma)
    return {"price": price, "rho_r": -B * price, "gamma_r": B * B * price,
            "duration": B, "convexity": B * B}


def cir_bond_option(r0, t_option, t_bond, strike, kappa, theta, sigma,
                    option_type=OptionType.CALL):
    """European option on a CIR zero-coupon bond (CIR 1985, exact).

    Option expires at ``t_option`` on a bond maturing at ``t_bond`` (``> t_option``),
    struck at ``strike`` on the bond price. Uses the noncentral chi-square
    formula: with ``g = sqrt(kappa^2 + 2 sigma^2)`` and the CIR affine
    ``A, B`` over ``t_bond - t_option``, the call is

        P(0,t_bond) X2(...; nc1) - strike P(0,t_option) X2(...; nc2),

    with critical rate ``r* = ln(A/strike)/B``. Puts follow from put-call parity.
    """
    ot = _coerce_type(option_type)
    if not (0 < t_option < t_bond):
        raise ValueError("require 0 < t_option < t_bond")
    if r0 < 0:
        raise ValueError("r0 must be non-negative in CIR")
    if kappa <= 0 or theta < 0 or sigma <= 0:
        raise ValueError("require kappa > 0, theta >= 0, sigma > 0")
    PS = cir_zero_coupon_bond(r0, t_option, kappa, theta, sigma)
    PT = cir_zero_coupon_bond(r0, t_bond, kappa, theta, sigma)

    g = math.sqrt(kappa * kappa + 2.0 * sigma * sigma)
    emgs = math.exp(g * t_option) - 1.0
    phi = 2.0 * g / (sigma * sigma * emgs)
    psi = (kappa + g) / (sigma * sigma)
    A_TS, B_TS = _AB(t_bond - t_option, kappa, theta, sigma)
    rstar = math.log(A_TS / strike) / B_TS
    df = 4.0 * kappa * theta / (sigma * sigma)
    common = 2.0 * phi * phi * r0 * math.exp(g * t_option)
    x1 = 2.0 * rstar * (phi + psi + B_TS)
    nc1 = common / (phi + psi + B_TS)
    x2 = 2.0 * rstar * (phi + psi)
    nc2 = common / (phi + psi)
    call = (PT * noncentral_chisq_cdf(x1, df, nc1)
            - strike * PS * noncentral_chisq_cdf(x2, df, nc2))
    if ot is OptionType.CALL:
        return call
    # Put-call parity: C - P = P(0,t_bond) - strike P(0,t_option).
    return call - PT + strike * PS


def cir_coupon_bond_option(r0, t_option, cashflows, strike, kappa, theta, sigma,
                           option_type=OptionType.CALL):
    """European option on a coupon bond under CIR (Jamshidian decomposition).

    ``cashflows`` is ``[(t_i, c_i), ...]`` with ``t_i > t_option``. The CIR bond
    is monotone decreasing in ``r0``, so Jamshidian applies: solve for the
    critical rate ``r*`` where the coupon bond's value at expiry equals
    ``strike``, then sum the ``c_i``-weighted CIR zero-coupon-bond options
    (:func:`cir_bond_option`) struck at ``K_i = P(t_option, t_i | r*)``. Exact.
    """
    ot = _coerce_type(option_type)
    cfs = sorted(cashflows)
    if not cfs or any(ti <= t_option for ti, _ in cfs):
        raise ValueError("all cashflow times must exceed t_option")
    if r0 < 0:
        raise ValueError("r0 must be non-negative in CIR")

    def bond_value_at(r):
        return sum(c * cir_zero_coupon_bond(r, ti - t_option, kappa, theta, sigma)
                   for ti, c in cfs)

    # CIR keeps r >= 0; bracket the critical rate on [0, hi].
    lo, hi = 0.0, 1.0
    while bond_value_at(hi) > strike and hi < 50.0:
        hi += 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if bond_value_at(mid) > strike:
            lo = mid
        else:
            hi = mid
    rstar = 0.5 * (lo + hi)

    total = 0.0
    for ti, c in cfs:
        Ki = cir_zero_coupon_bond(rstar, ti - t_option, kappa, theta, sigma)
        total += c * cir_bond_option(r0, t_option, ti, Ki, kappa, theta, sigma,
                                     ot)
    return total


def cir_swaption(r0, expiry, pay_times, fixed_rate, kappa, theta, sigma,
                 payer=True, notional=1.0):
    """European swaption under CIR via the coupon-bond-option identity (exact).

    A payer swaption is a put on the fixed-leg coupon bond struck at the
    notional; a receiver is a call. Priced by :func:`cir_coupon_bond_option`.
    ``pay_times`` are the fixed-leg payment dates (all ``> expiry``); accruals
    are the gaps, the first measured from ``expiry``.
    """
    times = sorted(pay_times)
    if not times or any(t <= expiry for t in times):
        raise ValueError("all pay_times must exceed expiry")
    prev = expiry
    cfs = []
    for i, ti in enumerate(times):
        tau = ti - prev
        prev = ti
        c = fixed_rate * tau * notional
        if i == len(times) - 1:
            c += notional
        cfs.append((ti, c))
    ot = OptionType.PUT if payer else OptionType.CALL
    return cir_coupon_bond_option(r0, expiry, cfs, notional, kappa, theta,
                                  sigma, ot)
