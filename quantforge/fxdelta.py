"""FX option delta-space quoting conventions.

FX vol markets quote by *delta*, not strike: a smile is given as an at-the-money
vol plus 25-delta (and 10-delta) risk reversals and butterflies. This module
converts between strike and delta and builds the market pillars:

  * ``strike_from_delta`` / ``delta_from_strike`` in the forward- or spot-delta
    convention, with an optional premium adjustment (for premium-in-foreign
    currency pairs);
  * ``atm_dns_strike`` -- the delta-neutral-straddle ATM strike
    ``K = F exp(0.5 sigma^2 t)`` (the market ATM convention);
  * ``rr_bf_to_pillars`` -- turn (ATM, RR, BF) at a delta into the call/put
    pillar vols and their strikes.

Pure standard library; uses the accurate inverse-normal from
:mod:`quantforge.mathfns`.
"""

import math

from .mathfns import norm_cdf, norm_ppf


def atm_dns_strike(F, t, sigma):
    """Delta-neutral-straddle ATM strike ``K = F exp(0.5 sigma^2 t)``.

    The strike at which a straddle has zero (forward) delta -- the market's
    standard ATM quote for most currency pairs.
    """
    return F * math.exp(0.5 * sigma * sigma * t)


def strike_from_delta(F, t, sigma, delta, is_call, premium_adjusted=False,
                      spot_delta=False, r_for=0.0):
    """Strike with the given delta.

    Args:
        delta: the target delta magnitude convention -- pass the signed delta
            (call > 0, put < 0), e.g. ``0.25`` for a 25-delta call, ``-0.25``
            for a 25-delta put.
        is_call: whether the option is a call.
        premium_adjusted: use the premium-adjusted delta convention (delta net
            of the option premium, standard for premium-in-foreign pairs).
        spot_delta: if True the delta is a spot delta (discounted by ``r_for``);
            otherwise a forward delta.

    Returns the strike ``K``.
    """
    vsqrt = sigma * math.sqrt(t)
    # Undo the foreign discounting for a spot delta to get the forward delta.
    fdelta = delta
    if spot_delta:
        fdelta = delta * math.exp(r_for * t)

    if not premium_adjusted:
        # forward delta_call = N(d1); delta_put = -N(-d1) = N(d1) - 1.
        if is_call:
            d1 = norm_ppf(fdelta)
        else:
            d1 = norm_ppf(1.0 + fdelta)   # put delta negative
        return F * math.exp(-d1 * vsqrt + 0.5 * vsqrt * vsqrt)

    # Premium-adjusted: delta = phi * (K/F) N(phi d2). Solve for d2 by Newton.
    phi = 1.0 if is_call else -1.0
    target = abs(fdelta)
    d2 = phi * norm_ppf(target) if is_call else -norm_ppf(target)
    for _ in range(100):
        KoverF = math.exp(-d2 * vsqrt - 0.5 * vsqrt * vsqrt)
        f = KoverF * norm_cdf(phi * d2) - target
        # derivative wrt d2.
        dKoverF = -vsqrt * KoverF
        from .mathfns import norm_pdf
        df = dKoverF * norm_cdf(phi * d2) + KoverF * phi * norm_pdf(phi * d2)
        if abs(df) < 1e-14:
            break
        step = f / df
        d2 -= step
        if abs(step) < 1e-12:
            break
    return F * math.exp(-d2 * vsqrt - 0.5 * vsqrt * vsqrt)


def delta_from_strike(F, t, sigma, K, is_call, premium_adjusted=False,
                      spot_delta=False, r_for=0.0):
    """Forward (or spot) delta of an option struck at ``K``.

    The inverse of :func:`strike_from_delta`; returns the signed delta.
    """
    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log(F / K) + 0.5 * vsqrt * vsqrt) / vsqrt
    d2 = d1 - vsqrt
    phi = 1.0 if is_call else -1.0
    if not premium_adjusted:
        delta = phi * norm_cdf(phi * d1)
    else:
        delta = phi * (K / F) * norm_cdf(phi * d2)
    if spot_delta:
        delta *= math.exp(-r_for * t)
    return delta


def rr_bf_to_pillars(F, t, atm, rr, bf, call_delta=0.25):
    """Convert (ATM, risk-reversal, butterfly) quotes to smile pillars.

    Returns ``(K_put, sigma_put, K_atm, sigma_atm, K_call, sigma_call)`` -- the
    three market pillar strikes and vols. Uses the standard smile-implied-from-
    quotes relations ``sigma_25c = atm + bf + rr/2``, ``sigma_25p = atm + bf -
    rr/2`` and the delta-neutral ATM strike.
    """
    sigma_call = atm + bf + 0.5 * rr
    sigma_put = atm + bf - 0.5 * rr
    K_atm = atm_dns_strike(F, t, atm)
    K_call = strike_from_delta(F, t, sigma_call, call_delta, True)
    K_put = strike_from_delta(F, t, sigma_put, -call_delta, False)
    return (K_put, sigma_put, K_atm, atm, K_call, sigma_call)
