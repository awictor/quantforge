"""Constant-Elasticity-of-Variance (CEV) option pricing.

The CEV model (Cox 1975) lets volatility depend on the spot level:

    dS = (r - q) S dt + delta * S^beta dW,   0 <= beta < 1 (equity leverage),

so the local vol ``delta * S^{beta-1}`` rises as the price falls (the leverage
effect / skew). Schroder (1989) gives a closed form for the European call in
terms of the **noncentral chi-square** distribution.

We implement the noncentral chi-square CDF from scratch (a Poisson-weighted sum
of central chi-square CDFs, i.e. regularized lower incomplete gammas), so the
whole model is dependency-free. The tests check the beta -> 1 limit against
Black-Scholes and put-call parity.
"""

import math

from .bsm import OptionType, _coerce_type, _validate


def _gammainc_lower_reg(a: float, x: float) -> float:
    """Regularized lower incomplete gamma P(a, x) = gamma(a, x) / Gamma(a).

    Series expansion for x < a+1, continued fraction otherwise (Numerical
    Recipes). Accurate to ~1e-12 over the range CEV needs.
    """
    if x < 0 or a <= 0:
        raise ValueError("require x >= 0 and a > 0")
    if x == 0:
        return 0.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        # Series representation.
        ap = a
        term = 1.0 / a
        total = term
        for _ in range(1000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return total * math.exp(-x + a * math.log(x) - gln)
    else:
        # Continued fraction (Lentz).
        tiny = 1e-300
        b = x + 1.0 - a
        c = 1.0 / tiny
        d = 1.0 / b
        h = d
        for i in range(1, 1000):
            an = -i * (i - a)
            b += 2.0
            d = an * d + b
            if abs(d) < tiny:
                d = tiny
            c = b + an / c
            if abs(c) < tiny:
                c = tiny
            d = 1.0 / d
            delta = d * c
            h *= delta
            if abs(delta - 1.0) < 1e-15:
                break
        q = math.exp(-x + a * math.log(x) - gln) * h
        return 1.0 - q


def noncentral_chisq_cdf(x: float, k: float, lam: float) -> float:
    """Noncentral chi-square CDF at ``x`` with ``k`` dof and noncentrality ``lam``.

    A Poisson(lam/2)-weighted sum of central chi-square CDFs:
    ``F(x; k, lam) = sum_j pois(j; lam/2) * P((k+2j)/2, x/2)``. The summation
    starts at the Poisson mode ``j0 = floor(lam/2)`` and expands outward, so it
    stays numerically stable even when ``lam`` is large (the naive j=0 start
    underflows because ``e^{-lam/2}`` is zero to machine precision).
    """
    if x <= 0:
        return 0.0
    if lam < 0 or k <= 0:
        raise ValueError("require k > 0 and lam >= 0")
    half_lam = 0.5 * lam
    half_x = 0.5 * x
    if half_lam == 0.0:
        return _gammainc_lower_reg(0.5 * k, half_x)

    j0 = int(half_lam)
    # log Poisson weight at the mode: -lam/2 + j0 ln(lam/2) - ln(j0!).
    log_w0 = -half_lam + j0 * math.log(half_lam) - math.lgamma(j0 + 1.0)
    w0 = math.exp(log_w0)

    total = w0 * _gammainc_lower_reg(0.5 * k + j0, half_x)

    # Expand upward from the mode.
    w = w0
    j = j0
    for _ in range(100000):
        j += 1
        w *= half_lam / j
        term = w * _gammainc_lower_reg(0.5 * k + j, half_x)
        total += term
        if w < 1e-18 and j > half_lam:
            break

    # Expand downward from the mode.
    w = w0
    j = j0
    while j > 0:
        w *= j / half_lam
        j -= 1
        total += w * _gammainc_lower_reg(0.5 * k + j, half_x)
        if w < 1e-18:
            break
    return min(max(total, 0.0), 1.0)


def cev_price(S, K, t, r, sigma, beta, option_type=OptionType.CALL, q=0.0):
    """European CEV option price (Schroder 1989), for 0 <= beta < 1.

    Args:
        sigma: the CEV volatility level, calibrated so that at ``S`` the
            instantaneous lognormal vol equals ``sigma`` — i.e. the scale
            ``delta = sigma * S^{1-beta}``. This makes ``beta`` control only the
            skew, with ``sigma`` comparable to a Black-Scholes vol.
        beta: elasticity in [0, 1). beta = 1 recovers Black-Scholes (handled by
            a limit); beta = 0 is the Bachelier-like absolute-diffusion case.
        q: continuous dividend yield.

    Puts are obtained by put-call parity.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if not (0.0 <= beta < 1.0):
        raise ValueError("beta must be in [0, 1)")
    if t == 0 or sigma == 0:
        fwd = S * math.exp((r - q) * t)
        disc = math.exp(-r * t)
        payoff = max(fwd - K, 0.0) if ot is OptionType.CALL else max(K - fwd, 0.0)
        return disc * payoff

    mu = r - q
    # Scale so the ATM instantaneous vol matches sigma: delta = sigma S^{1-beta}.
    delta = sigma * S ** (1.0 - beta)
    one_mb = 1.0 - beta

    # Schroder / Hull noncentral chi-square parameters.
    if abs(mu) > 1e-12:
        v = delta * delta / (2.0 * mu * (beta - 1.0)) * (
            math.exp(2.0 * mu * (beta - 1.0) * t) - 1.0)
    else:
        v = delta * delta * t   # mu -> 0 limit
    a = (K * math.exp(-mu * t)) ** (2.0 * one_mb) / (one_mb * one_mb * v)
    bb = 1.0 / one_mb
    c = S ** (2.0 * one_mb) / (one_mb * one_mb * v)

    disc = math.exp(-r * t)
    call = (S * math.exp(-q * t) * (1.0 - noncentral_chisq_cdf(a, bb + 2.0, c))
            - K * disc * noncentral_chisq_cdf(c, bb, a))

    if ot is OptionType.CALL:
        return call
    # Put-call parity: P = C - S e^{-qt} + K e^{-rt}.
    return call - S * math.exp(-q * t) + K * disc
