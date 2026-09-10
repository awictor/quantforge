"""Smile-aware delta adjustments: sticky-strike vs sticky-delta (sticky-moneyness).

The Black-Scholes delta assumes volatility is fixed as the spot moves. In a
real market the implied-vol smile moves with the spot, so the *effective*
delta of an option includes a second term, ``vega * d(sigma)/d(spot)`` — how
the option's own implied vol shifts as the underlying moves. Two standard
regimes describe that shift:

  * **sticky-strike**: the vol at each fixed strike stays put as spot moves, so
    ``d(sigma_K)/d(spot) = 0`` and the effective delta equals the BS delta.
    The smile still contributes through the skew when you reprice, but the
    instantaneous hedge ratio is unchanged.
  * **sticky-delta** (sticky-moneyness): the smile is a function of moneyness
    (e.g. ``k = ln(K/F)``), so it rides along with the spot. The vol at a fixed
    strike then changes as ``d(sigma_K)/d(spot) = -(d(sigma)/dk) / S`` and the
    effective delta picks up a ``vega`` correction.

``smile_delta`` returns the adjusted delta given a local skew slope
``dsigma/dk``; ``skew_slope`` estimates that slope by finite-differencing a
vol-smile function you supply.
"""

import math
from enum import Enum

from .bsm import delta as bs_delta, vega as bs_vega, OptionType, _coerce_type, _validate


class StickyRule(str, Enum):
    STRIKE = "strike"
    DELTA = "delta"        # a.k.a. sticky-moneyness


def skew_slope(smile_fn, K, F, h=None):
    """Estimate d(sigma)/dk at strike ``K`` via central finite difference.

    Args:
        smile_fn: callable ``sigma(K)`` returning implied vol for a strike.
        K: strike at which to measure the slope.
        F: forward (used to convert to log-moneyness k = ln(K/F)).
        h: bump in ``k`` space; defaults to a small fraction.

    Returns d(sigma)/dk where k = ln(K/F).
    """
    if h is None:
        h = 1e-4
    # k = ln(K/F)  =>  K = F * exp(k). Bump k by +/- h.
    k = math.log(K / F)
    K_up = F * math.exp(k + h)
    K_dn = F * math.exp(k - h)
    return (smile_fn(K_up) - smile_fn(K_dn)) / (2.0 * h)


def smile_delta(S, K, t, r, sigma, dsigma_dk=0.0,
                option_type=OptionType.CALL, b=None,
                sticky=StickyRule.DELTA) -> float:
    """Effective (smile-adjusted) delta of an option.

    Args:
        sigma: the option's current implied volatility.
        dsigma_dk: local skew slope d(sigma)/dk at this strike, where
            k = ln(K/F). Only used under the sticky-delta rule.
        sticky: STRIKE (delta == BS delta) or DELTA (add the vega/skew term).

    Under sticky-delta the smile is a function of moneyness, so a 1-unit rise in
    spot lowers the log-moneyness of a fixed strike by 1/S, shifting its vol by
    ``-(dsigma/dk)/S``. The effective delta is therefore

        delta_eff = delta_BS + vega * d(sigma)/d(spot)
                  = delta_BS - vega * (dsigma/dk) / S.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    rule = StickyRule(sticky)
    d_bs = bs_delta(S, K, t, r, sigma, ot, b)
    if rule is StickyRule.STRIKE:
        return d_bs
    v = bs_vega(S, K, t, r, sigma, b)
    return d_bs - v * dsigma_dk / S


def smile_delta_from_smile(S, K, t, r, smile_fn, F=None,
                           option_type=OptionType.CALL, b=None) -> float:
    """Sticky-delta effective delta computed directly from a smile function.

    Reads the option's vol as ``smile_fn(K)``, estimates the local skew slope by
    finite difference, and returns the adjusted delta. ``F`` defaults to the
    carry-implied forward ``S * exp(b * t)``.
    """
    if b is None:
        b = r
    if F is None:
        F = S * math.exp(b * t)
    sigma = smile_fn(K)
    slope = skew_slope(smile_fn, K, F)
    return smile_delta(S, K, t, r, sigma, dsigma_dk=slope,
                       option_type=option_type, b=b, sticky=StickyRule.DELTA)
