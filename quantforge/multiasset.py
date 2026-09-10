"""Two-asset options: exchange (Margrabe), spread (Kirk), and basket (moment-match).

These price options whose payoff depends on two correlated lognormal assets:

  * **exchange** — the right to swap asset 2 for asset 1, payoff
    ``max(S1 - S2, 0)``. Margrabe (1978) gives an exact closed form: a
    Black-Scholes call with zero strike and an effective volatility
    ``sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)``.
  * **spread** — payoff ``max(S1 - S2 - K, 0)`` for a non-zero strike. No exact
    closed form exists; Kirk's (1995) approximation reduces it to a Margrabe-
    style formula and is accurate for moderate strikes.
  * **basket** — a call on ``w1 S1 + w2 S2``. The sum of lognormals is not
    lognormal, so we moment-match it to a single lognormal (Levy) and apply
    Black-Scholes.

All take continuous dividend yields ``q1, q2`` and the correlation ``rho``.
"""

import math

from .mathfns import norm_cdf
from .bsm import price as bsm_price, OptionType, _coerce_type


def exchange_option(S1, S2, t, sigma1, sigma2, rho, q1=0.0, q2=0.0) -> float:
    """Margrabe option to exchange asset 2 for asset 1: payoff max(S1 - S2, 0).

    Exact closed form; independent of the risk-free rate (the two assets'
    financing cancels), depending only on the dividend yields.
    """
    if S1 <= 0 or S2 <= 0:
        raise ValueError("prices must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    sigma = math.sqrt(sigma1 * sigma1 - 2.0 * rho * sigma1 * sigma2 + sigma2 * sigma2)
    if t == 0 or sigma == 0:
        return max(S1 * math.exp(-q1 * t) - S2 * math.exp(-q2 * t), 0.0)
    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log((S1 * math.exp(-q1 * t)) / (S2 * math.exp(-q2 * t)))
          + 0.5 * sigma * sigma * t) / vsqrt
    d2 = d1 - vsqrt
    return S1 * math.exp(-q1 * t) * norm_cdf(d1) - S2 * math.exp(-q2 * t) * norm_cdf(d2)


def spread_option(S1, S2, K, t, r, sigma1, sigma2, rho,
                  q1=0.0, q2=0.0, option_type=OptionType.CALL) -> float:
    """Kirk (1995) approximation for a spread option: payoff max(S1 - S2 - K, 0).

    Reduces to an exact Margrabe formula when K = 0. Puts follow from parity on
    the spread ``S1 - S2``.
    """
    ot = _coerce_type(option_type)
    if S1 <= 0 or S2 <= 0:
        raise ValueError("prices must be positive")

    # Forwards of each asset.
    F1 = S1 * math.exp((r - q1) * t)
    F2 = S2 * math.exp((r - q2) * t)
    disc = math.exp(-r * t)

    if t == 0:
        payoff = (max(S1 - S2 - K, 0.0) if ot is OptionType.CALL
                  else max(K - (S1 - S2), 0.0))
        return payoff

    # Kirk: treat (F2 + K) as a single lognormal asset with a blended vol.
    denom = F2 + K
    ratio = F2 / denom
    sigma_k = math.sqrt(sigma1 * sigma1
                        - 2.0 * rho * sigma1 * sigma2 * ratio
                        + sigma2 * sigma2 * ratio * ratio)
    vsqrt = sigma_k * math.sqrt(t)
    d1 = (math.log(F1 / denom) + 0.5 * sigma_k * sigma_k * t) / vsqrt
    d2 = d1 - vsqrt
    call = disc * (F1 * norm_cdf(d1) - denom * norm_cdf(d2))
    if ot is OptionType.CALL:
        return call
    # Put via parity: C - P = disc*(F1 - F2 - K).
    return call - disc * (F1 - F2 - K)


def basket_option(spots, weights, K, t, r, sigmas, corr, q=None,
                  option_type=OptionType.CALL) -> float:
    """Two-asset basket call/put on ``w1 S1 + w2 S2`` via lognormal moment match.

    Args:
        spots: (S1, S2). weights: (w1, w2). sigmas: (sigma1, sigma2).
        corr: correlation between the two assets.
        q: optional (q1, q2) dividend yields; defaults to zeros.

    Matches the basket forward's first two moments to a single lognormal (Levy)
    and prices with Black-Scholes. Exact for a single asset; an approximation
    for the sum.
    """
    ot = _coerce_type(option_type)
    if len(spots) != 2 or len(weights) != 2 or len(sigmas) != 2:
        raise ValueError("basket_option handles exactly two assets")
    if q is None:
        q = (0.0, 0.0)
    S1, S2 = spots
    w1, w2 = weights
    s1, s2 = sigmas
    q1, q2 = q

    F1 = w1 * S1 * math.exp((r - q1) * t)
    F2 = w2 * S2 * math.exp((r - q2) * t)
    M1 = F1 + F2  # first moment of the basket forward

    if M1 <= 0:
        raise ValueError("basket forward must be positive")
    if t == 0:
        basket = w1 * S1 + w2 * S2
        return max(basket - K, 0.0) if ot is OptionType.CALL else max(K - basket, 0.0)

    # Second moment of the basket forward (cross term carries the correlation).
    M2 = (F1 * F1 * math.exp(s1 * s1 * t)
          + 2.0 * F1 * F2 * math.exp(corr * s1 * s2 * t)
          + F2 * F2 * math.exp(s2 * s2 * t))
    # Implied lognormal vol from the two moments.
    sigma_b = math.sqrt(math.log(M2 / (M1 * M1)) / t)

    # Black-Scholes on the basket forward, discounted.
    disc = math.exp(-r * t)
    vsqrt = sigma_b * math.sqrt(t)
    d1 = (math.log(M1 / K) + 0.5 * sigma_b * sigma_b * t) / vsqrt
    d2 = d1 - vsqrt
    if ot is OptionType.CALL:
        return disc * (M1 * norm_cdf(d1) - K * norm_cdf(d2))
    return disc * (K * norm_cdf(-d2) - M1 * norm_cdf(-d1))
