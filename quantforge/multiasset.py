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
from .american import _bivariate_normal


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


def _rainbow_mc(S1, S2, K, t, r, sigma1, sigma2, rho, q1, q2, kind, ot,
                n_paths, antithetic, seed):
    """Shared correlated-GBM Monte Carlo for best-of / worst-of payoffs."""
    import random
    rng = random.Random(seed)
    d1 = (r - q1 - 0.5 * sigma1 * sigma1) * t
    d2 = (r - q2 - 0.5 * sigma2 * sigma2) * t
    v1 = sigma1 * math.sqrt(t)
    v2 = sigma2 * math.sqrt(t)
    corr2 = math.sqrt(1.0 - rho * rho)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    def payoff(z1, z2):
        a = S1 * math.exp(d1 + v1 * z1)
        b_ = S2 * math.exp(d2 + v2 * (rho * z1 + corr2 * z2))
        chosen = max(a, b_) if kind == "best" else min(a, b_)
        return max(sign * (chosen - K), 0.0)

    total = 0.0
    n = n_paths // 2 if antithetic else n_paths
    count = 0
    for _ in range(n):
        z1 = rng.gauss(0.0, 1.0)
        z2 = rng.gauss(0.0, 1.0)
        total += payoff(z1, z2)
        count += 1
        if antithetic:
            total += payoff(-z1, -z2)
            count += 1
    return disc * total / count


def _stulz_min_call(S1, S2, K, t, r, sigma1, sigma2, rho, q1, q2):
    """Exact Stulz (1982) price of a call on the minimum of two assets.

    ``max(min(S1, S2) - K, 0)``. With
    ``sigma = sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)`` the spread vol,
    ``rho1 = (sigma1 - rho sigma2)/sigma``, ``rho2 = (sigma2 - rho sigma1)/sigma``,

        d  = (ln(S1/S2) + (q2 - q1 + sigma^2/2) t) / (sigma sqrt(t))
        yi = (ln(Si/K) + (r - qi + sigma_i^2/2) t) / (sigma_i sqrt(t))

    the price is (Stulz eq. for the min):

        S1 e^{-q1 t} M(y1, -d; -rho1)
      + S2 e^{-q2 t} M(y2, d - sigma sqrt(t); -rho2)
      - K e^{-r t} M(y1 - sigma1 sqrt(t), y2 - sigma2 sqrt(t); rho)

    where ``M`` is the standardized bivariate-normal CDF.
    """
    st = math.sqrt(t)
    sig = math.sqrt(sigma1 * sigma1 - 2.0 * rho * sigma1 * sigma2 + sigma2 * sigma2)
    disc = math.exp(-r * t)
    if sig < 1e-12:
        # Perfectly correlated with equal vols: min is a single lognormal.
        c1 = bsm_price(S1, K, t, r, sigma1, OptionType.CALL, b=r - q1)
        c2 = bsm_price(S2, K, t, r, sigma2, OptionType.CALL, b=r - q2)
        return min(c1, c2)
    rho1 = (sigma1 - rho * sigma2) / sig
    rho2 = (sigma2 - rho * sigma1) / sig
    d = (math.log(S1 / S2) + (q2 - q1 + 0.5 * sig * sig) * t) / (sig * st)
    y1 = (math.log(S1 / K) + (r - q1 + 0.5 * sigma1 * sigma1) * t) / (sigma1 * st)
    y2 = (math.log(S2 / K) + (r - q2 + 0.5 * sigma2 * sigma2) * t) / (sigma2 * st)
    return (S1 * math.exp(-q1 * t) * _bivariate_normal(y1, -d, -rho1)
            + S2 * math.exp(-q2 * t) * _bivariate_normal(y2, d - sig * st, -rho2)
            - K * disc * _bivariate_normal(y1 - sigma1 * st, y2 - sigma2 * st, rho))


def best_of_call_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0):
    """Exact Stulz (1982) price of a call on the maximum of two assets.

    ``max(max(S1, S2) - K, 0)``. Uses the Stulz identity
    ``C_max + C_min = c(S1) + c(S2)`` (both vanilla calls at strike ``K``), so
    ``C_max = c(S1) + c(S2) - C_min`` with the exact :func:`_stulz_min_call`.
    This is the closed-form cross-check for the Monte Carlo :func:`best_of_call`.
    """
    if S1 <= 0 or S2 <= 0 or K <= 0:
        raise ValueError("prices and strike must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")
    c1 = bsm_price(S1, K, t, r, sigma1, OptionType.CALL, b=r - q1)
    c2 = bsm_price(S2, K, t, r, sigma2, OptionType.CALL, b=r - q2)
    return c1 + c2 - _stulz_min_call(S1, S2, K, t, r, sigma1, sigma2, rho, q1, q2)


def worst_of_call_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0):
    """Exact Stulz (1982) price of a call on the minimum of two assets.

    ``max(min(S1, S2) - K, 0)``. Closed-form cross-check for the Monte Carlo
    :func:`worst_of_call`.
    """
    if S1 <= 0 or S2 <= 0 or K <= 0:
        raise ValueError("prices and strike must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")
    return _stulz_min_call(S1, S2, K, t, r, sigma1, sigma2, rho, q1, q2)


def two_asset_digital(S1, S2, K1, K2, t, r, sigma1, sigma2, rho,
                      cond1="above", cond2="above", q1=0.0, q2=0.0, cash=1.0):
    """Cash-or-nothing digital on two correlated assets (exact closed form).

    Pays ``cash`` at expiry iff both single-asset conditions hold: asset 1 is
    ``above`` (``S1_T > K1``) or ``below`` (``S1_T < K1``) its strike, and
    likewise for asset 2. Under the risk-neutral bivariate lognormal the price is

        cash * e^{-r t} * M(s1 d1, s2 d2; s1 s2 rho)

    where ``di = (ln(Si/Ki) + (r - qi - sigma_i^2/2) t) / (sigma_i sqrt(t))`` is
    the usual ``d2``, ``si = +1`` for an ``above`` condition and ``-1`` for a
    ``below`` one, and ``M`` is the standardized bivariate-normal CDF. Flipping a
    condition flips the sign of that ``d`` and of the correlation. The four
    quadrant prices sum to ``cash e^{-r t}`` (the conditions are exhaustive).

    Cross-checks a correlated-GBM Monte Carlo.
    """
    if S1 <= 0 or S2 <= 0 or K1 <= 0 or K2 <= 0:
        raise ValueError("prices and strikes must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")
    c1 = str(cond1).lower()
    c2 = str(cond2).lower()
    if c1 not in ("above", "below") or c2 not in ("above", "below"):
        raise ValueError("cond1/cond2 must be 'above' or 'below'")
    st = math.sqrt(t)
    d1 = (math.log(S1 / K1) + (r - q1 - 0.5 * sigma1 * sigma1) * t) / (sigma1 * st)
    d2 = (math.log(S2 / K2) + (r - q2 - 0.5 * sigma2 * sigma2) * t) / (sigma2 * st)
    s1 = 1.0 if c1 == "above" else -1.0
    s2 = 1.0 if c2 == "above" else -1.0
    disc = math.exp(-r * t)
    return cash * disc * _bivariate_normal(s1 * d1, s2 * d2, s1 * s2 * rho)


def two_asset_asset_or_nothing(S1, S2, K1, K2, t, r, sigma1, sigma2, rho,
                               cond1="above", cond2="above", q1=0.0, q2=0.0):
    """Asset-or-nothing digital paying ``S1_T`` iff both conditions hold.

    Pays the *first* asset's terminal value at expiry iff asset 1 is
    ``above``/``below`` ``K1`` and asset 2 is ``above``/``below`` ``K2``. Pricing
    under the asset-1 (share) measure -- where ``S1`` is the numeraire and asset
    1's drift gains ``sigma1^2`` while asset 2's shock inherits an extra
    ``rho sigma1 sqrt(t)`` -- gives

        S1 e^{-q1 t} * M(s1 a1, s2 a2; s1 s2 rho),

        a1 = (ln(S1/K1) + (r - q1 + sigma1^2/2) t) / (sigma1 sqrt t)
        a2 = (ln(S2/K2) + (r - q2 - sigma2^2/2) t) / (sigma2 sqrt t)
             + rho sigma1 sqrt(t)

    with ``si = +1`` for an ``above`` condition, ``-1`` for ``below``. The four
    quadrant prices sum to the discounted forward ``S1 e^{-q1 t}`` (asset 1 is
    always delivered, on some quadrant). Cross-checks a correlated-GBM Monte
    Carlo. To pay asset 2 instead, swap the two assets in the call.
    """
    if S1 <= 0 or S2 <= 0 or K1 <= 0 or K2 <= 0:
        raise ValueError("prices and strikes must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")
    c1 = str(cond1).lower()
    c2 = str(cond2).lower()
    if c1 not in ("above", "below") or c2 not in ("above", "below"):
        raise ValueError("cond1/cond2 must be 'above' or 'below'")
    st = math.sqrt(t)
    a1 = (math.log(S1 / K1) + (r - q1 + 0.5 * sigma1 * sigma1) * t) / (sigma1 * st)
    a2 = (math.log(S2 / K2) + (r - q2 - 0.5 * sigma2 * sigma2) * t) / (sigma2 * st) \
        + rho * sigma1 * st
    s1 = 1.0 if c1 == "above" else -1.0
    s2 = 1.0 if c2 == "above" else -1.0
    return S1 * math.exp(-q1 * t) * _bivariate_normal(s1 * a1, s2 * a2, s1 * s2 * rho)


def _disc_expected_min(S1, S2, t, r, sigma1, sigma2, rho, q1, q2):
    """Discounted risk-neutral expectation of ``min(S1_T, S2_T)``.

    ``min(a, b) = b - max(b - a, 0)``, so the discounted expectation of the min
    is the discounted forward of asset 2 minus the Margrabe value of the option
    to exchange asset 1 for asset 2 (payoff ``max(S2 - S1, 0)``):

        disc E[min] = S2 e^{-q2 t} - exchange_option(S2, S1; ...).
    """
    return S2 * math.exp(-q2 * t) - exchange_option(S2, S1, t, sigma2, sigma1,
                                                     rho, q2, q1)


def worst_of_put_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0):
    """Exact price of a put on the minimum of two assets: ``max(K - min(S1,S2), 0)``.

    By put-call parity on the rainbow, a put and call on the same underlying
    (here ``min(S1, S2)``) satisfy ``C - P = disc E[min] - K e^{-r t}``, so

        P_min = C_min - disc E[min] + K e^{-r t}

    with the exact :func:`worst_of_call_closed` and :func:`_disc_expected_min`.
    Closed-form cross-check for the Monte Carlo :func:`worst_of_call` put.
    """
    if S1 <= 0 or S2 <= 0 or K <= 0:
        raise ValueError("prices and strike must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")
    cmin = worst_of_call_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1, q2)
    demin = _disc_expected_min(S1, S2, t, r, sigma1, sigma2, rho, q1, q2)
    return cmin - demin + K * math.exp(-r * t)


def best_of_put_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0):
    """Exact price of a put on the maximum of two assets: ``max(K - max(S1,S2), 0)``.

        P_max = C_max - disc E[max] + K e^{-r t},

    where ``disc E[max] = S1 e^{-q1 t} + S2 e^{-q2 t} - disc E[min]`` (the two
    forwards less the discounted expected min). Uses the exact
    :func:`best_of_call_closed`. Cross-check for the Monte Carlo
    :func:`best_of_call` put.
    """
    if S1 <= 0 or S2 <= 0 or K <= 0:
        raise ValueError("prices and strike must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")
    cmax = best_of_call_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1, q2)
    demin = _disc_expected_min(S1, S2, t, r, sigma1, sigma2, rho, q1, q2)
    demax = S1 * math.exp(-q1 * t) + S2 * math.exp(-q2 * t) - demin
    return cmax - demax + K * math.exp(-r * t)


def best_of_call(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0,
                 option_type=OptionType.CALL, n_paths=100_000, antithetic=True,
                 seed=None):
    """Option on the maximum of two assets: payoff max(max(S1,S2) - K, 0) (call).

    Monte Carlo on correlated GBM. Best-of and worst-of calls satisfy
    ``best + worst = call(S1) + call(S2)`` at the same strike (Stulz), which the
    tests check.
    """
    ot = _coerce_type(option_type)
    return _rainbow_mc(S1, S2, K, t, r, sigma1, sigma2, rho, q1, q2, "best", ot,
                       n_paths, antithetic, seed)


def worst_of_call(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0,
                  option_type=OptionType.CALL, n_paths=100_000, antithetic=True,
                  seed=None):
    """Option on the minimum of two assets: payoff max(min(S1,S2) - K, 0) (call)."""
    ot = _coerce_type(option_type)
    return _rainbow_mc(S1, S2, K, t, r, sigma1, sigma2, rho, q1, q2, "worst", ot,
                       n_paths, antithetic, seed)


def exchange_greeks(S1, S2, t, sigma1, sigma2, rho, q1=0.0, q2=0.0):
    """Greeks of a Margrabe exchange option (payoff max(S1 - S2, 0)) by FD.

    Returns a dict with the two spot deltas (``delta1`` = dV/dS1,
    ``delta2`` = dV/dS2), the two own-gammas (``gamma1``, ``gamma2``), the
    cross-gamma (``cross`` = d2V/dS1 dS2), and the correlation sensitivity
    (``corr_vega`` = dV/drho). All by central finite differences on the exact
    Margrabe formula.
    """
    def px(a=S1, b_=S2, s1=sigma1, s2=sigma2, rr=rho):
        return exchange_option(a, b_, t, s1, s2, rr, q1, q2)

    base = px()
    h1 = 1e-3 * S1
    h2 = 1e-3 * S2
    d1u, d1d = px(a=S1 + h1), px(a=S1 - h1)
    d2u, d2d = px(b_=S2 + h2), px(b_=S2 - h2)
    delta1 = (d1u - d1d) / (2.0 * h1)
    delta2 = (d2u - d2d) / (2.0 * h2)
    gamma1 = (d1u - 2.0 * base + d1d) / (h1 * h1)
    gamma2 = (d2u - 2.0 * base + d2d) / (h2 * h2)
    # Cross-gamma via the mixed central difference.
    pp = px(a=S1 + h1, b_=S2 + h2)
    pm = px(a=S1 + h1, b_=S2 - h2)
    mp = px(a=S1 - h1, b_=S2 + h2)
    mm = px(a=S1 - h1, b_=S2 - h2)
    cross = (pp - pm - mp + mm) / (4.0 * h1 * h2)

    hr = 1e-5
    corr_vega = (px(rr=min(rho + hr, 1.0 - 1e-9))
                 - px(rr=max(rho - hr, -1.0 + 1e-9))) / (2.0 * hr)

    return {"price": base, "delta1": delta1, "delta2": delta2,
            "gamma1": gamma1, "gamma2": gamma2, "cross": cross,
            "corr_vega": corr_vega}


def spread_greeks(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0,
                  option_type=OptionType.CALL):
    """Greeks of a Kirk spread option (payoff max(S1 - S2 - K, 0)) by FD.

    Returns a dict with the two spot deltas, own-gammas, the cross-gamma
    (``d2V/dS1 dS2``), and the correlation sensitivity (``corr_vega``). All by
    central finite differences on the Kirk approximation.
    """
    ot = _coerce_type(option_type)

    def px(a=S1, b_=S2, rr=rho):
        return spread_option(a, b_, K, t, r, sigma1, sigma2, rr, q1, q2, ot)

    base = px()
    h1 = 1e-3 * S1
    h2 = 1e-3 * S2
    d1u, d1d = px(a=S1 + h1), px(a=S1 - h1)
    d2u, d2d = px(b_=S2 + h2), px(b_=S2 - h2)
    delta1 = (d1u - d1d) / (2.0 * h1)
    delta2 = (d2u - d2d) / (2.0 * h2)
    gamma1 = (d1u - 2.0 * base + d1d) / (h1 * h1)
    gamma2 = (d2u - 2.0 * base + d2d) / (h2 * h2)
    pp = px(a=S1 + h1, b_=S2 + h2)
    pm = px(a=S1 + h1, b_=S2 - h2)
    mp = px(a=S1 - h1, b_=S2 + h2)
    mm = px(a=S1 - h1, b_=S2 - h2)
    cross = (pp - pm - mp + mm) / (4.0 * h1 * h2)
    hr = 1e-5
    corr_vega = (px(rr=min(rho + hr, 1.0 - 1e-9))
                 - px(rr=max(rho - hr, -1.0 + 1e-9))) / (2.0 * hr)
    return {"price": base, "delta1": delta1, "delta2": delta2,
            "gamma1": gamma1, "gamma2": gamma2, "cross": cross,
            "corr_vega": corr_vega}


def basket_greeks(spots, weights, K, t, r, sigmas, corr, q=None,
                  option_type=OptionType.CALL):
    """Greeks of a two-asset basket option (Levy moment-match) by FD.

    Returns a dict with the two spot deltas, own-gammas, the cross-gamma, and
    the correlation sensitivity, all by central finite differences on
    :func:`basket_option`.
    """
    ot = _coerce_type(option_type)
    S1, S2 = spots

    def px(a=S1, b_=S2, rr=corr):
        return basket_option((a, b_), weights, K, t, r, sigmas, rr, q, ot)

    base = px()
    h1 = 1e-3 * S1
    h2 = 1e-3 * S2
    d1u, d1d = px(a=S1 + h1), px(a=S1 - h1)
    d2u, d2d = px(b_=S2 + h2), px(b_=S2 - h2)
    delta1 = (d1u - d1d) / (2.0 * h1)
    delta2 = (d2u - d2d) / (2.0 * h2)
    gamma1 = (d1u - 2.0 * base + d1d) / (h1 * h1)
    gamma2 = (d2u - 2.0 * base + d2d) / (h2 * h2)
    pp = px(a=S1 + h1, b_=S2 + h2)
    pm = px(a=S1 + h1, b_=S2 - h2)
    mp = px(a=S1 - h1, b_=S2 + h2)
    mm = px(a=S1 - h1, b_=S2 - h2)
    cross = (pp - pm - mp + mm) / (4.0 * h1 * h2)
    hr = 1e-5
    corr_vega = (px(rr=min(corr + hr, 1.0 - 1e-9))
                 - px(rr=max(corr - hr, -1.0 + 1e-9))) / (2.0 * hr)
    return {"price": base, "delta1": delta1, "delta2": delta2,
            "gamma1": gamma1, "gamma2": gamma2, "cross": cross,
            "corr_vega": corr_vega}


def implied_spread_correlation(target_price, S1, S2, K, t, r, sigma1, sigma2,
                               q1=0.0, q2=0.0, option_type=OptionType.CALL,
                               tol=1e-8, max_iter=100):
    """Back out the correlation implied by a spread-option market price (Kirk).

    The Kirk spread price is monotone decreasing in ``rho`` (higher correlation
    lowers the spread volatility), so a bisection on ``rho in (-1, 1)`` recovers
    the correlation consistent with the quote. Raises if the quote lies outside
    the price range spanned by ``rho = -1 .. 1``.
    """
    ot = _coerce_type(option_type)

    def px(rho):
        return spread_option(S1, S2, K, t, r, sigma1, sigma2, rho, q1, q2, ot)

    lo, hi = -0.999999, 0.999999
    p_lo, p_hi = px(lo), px(hi)      # p_lo is the highest price (rho=-1)
    if not (min(p_lo, p_hi) - 1e-10 <= target_price <= max(p_lo, p_hi) + 1e-10):
        raise ValueError(
            f"price {target_price} outside the rho-range [{p_hi:.6g}, {p_lo:.6g}]"
        )
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        pm = px(mid)
        if abs(pm - target_price) < tol:
            return mid
        # Price decreases in rho: if model price too high, raise rho.
        if pm > target_price:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
