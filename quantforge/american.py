"""Bjerksund-Stensland (2002) closed-form American option approximation.

The binomial tree in :mod:`quantforge.binomial` is exact in the limit but
O(steps^2). Bjerksund-Stensland is a closed-form approximation that prices an
American option in constant time by modeling early exercise with a flat
exercise boundary split into two time regions. It is accurate to a few basis
points against the tree for typical inputs and is the standard fast American
pricer on trading desks.

We implement the call formula directly; American puts use the Bjerksund-
Stensland put-call transformation

    P(S, K, r, b, sigma) = C(K, S, r - b, -b, sigma)

which is exact for American options. When early exercise has no value (an
American call with cost of carry ``b >= r``, i.e. no dividends), the price
collapses to the European Black-Scholes value, which we return directly.
"""

import math

from .mathfns import norm_cdf
from .bsm import OptionType, _coerce_type, _validate, price as bsm_price


def _phi(S, t, gamma, h, I, r, b, sigma):
    """Bjerksund-Stensland phi helper (Haug 2002 form)."""
    vsqrt = sigma * math.sqrt(t)
    lam = (-r + gamma * b + 0.5 * gamma * (gamma - 1.0) * sigma * sigma) * t
    kappa = 2.0 * b / (sigma * sigma) + (2.0 * gamma - 1.0)
    d = -(math.log(S / h) + (b + (gamma - 0.5) * sigma * sigma) * t) / vsqrt
    return (math.exp(lam) * (S ** gamma)
            * (norm_cdf(d)
               - ((I / S) ** kappa) * norm_cdf(d - 2.0 * math.log(I / S) / vsqrt)))


def _psi(S, t2, gamma, h, I2, I1, r, b, sigma, t1):
    """Bjerksund-Stensland psi helper (bivariate-normal term, Haug 2002 form).

    Signature matches Haug: psi(S, T2, gamma, h, I2, I1, t1).
    """
    cbnd = _bivariate_normal
    v1 = sigma * math.sqrt(t1)
    v2 = sigma * math.sqrt(t2)
    bt = b + (gamma - 0.5) * sigma * sigma
    kappa = 2.0 * b / (sigma * sigma) + (2.0 * gamma - 1.0)
    lam = -r + gamma * b + 0.5 * gamma * (gamma - 1.0) * sigma * sigma

    e1 = (math.log(S / I1) + bt * t1) / v1
    e2 = (math.log(I2 * I2 / (S * I1)) + bt * t1) / v1
    e3 = (math.log(S / I1) - bt * t1) / v1
    e4 = (math.log(I2 * I2 / (S * I1)) - bt * t1) / v1

    f1 = (math.log(S / h) + bt * t2) / v2
    f2 = (math.log(I2 * I2 / (S * h)) + bt * t2) / v2
    f3 = (math.log(I1 * I1 / (S * h)) + bt * t2) / v2
    f4 = (math.log(S * I1 * I1 / (h * I2 * I2)) + bt * t2) / v2

    rho = math.sqrt(t1 / t2)
    term = (S ** gamma) * math.exp(lam * t2)
    return term * (
        cbnd(-e1, -f1, rho)
        - ((I2 / S) ** kappa) * cbnd(-e2, -f2, rho)
        - ((I1 / S) ** kappa) * cbnd(-e3, -f3, -rho)
        + ((I1 / I2) ** kappa) * cbnd(-e4, -f4, -rho)
    )


def _bivariate_normal(a, b, rho):
    """Standardized bivariate normal CDF P(X<=a, Y<=b; corr=rho).

    Drezner-Wesolowsky single-integral form evaluated with a symmetric
    Gauss-Legendre rule. The tabulated weights below are pre-scaled by
    1/(2*pi), so the arc integral (substitution theta = asin(rho)*(t+1)/2,
    dtheta = asin(rho)/2 dt) contributes a factor of ``asin(rho)/2`` rather
    than the ``asin(rho)/(4*pi)`` used with unscaled weights.

    Accurate to ~1e-7 for |rho| up to ~0.98, which covers every correlation
    the Bjerksund-Stensland model produces (rho = sqrt(t1/t) ~ 0.79).
    """
    if a <= -1e10 or b <= -1e10:
        return 0.0
    if a >= 1e10:
        return norm_cdf(b)
    if b >= 1e10:
        return norm_cdf(a)
    if abs(rho) < 1e-12:
        return norm_cdf(a) * norm_cdf(b)

    # Symmetric Gauss-Legendre nodes on [0,1] with weights pre-scaled by 1/(2pi).
    x = [0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992]
    w = [0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042]

    h1, hk = a, b
    h12 = 0.5 * (h1 * h1 + hk * hk)
    hs = h1 * hk
    asr = math.asin(rho)
    bvn = 0.0
    for i in range(5):
        for sign in (-1.0, 1.0):
            sn = math.sin(asr * (sign * x[i] + 1.0) / 2.0)
            bvn += w[i] * math.exp((sn * hs - h12) / (1.0 - sn * sn))
    bvn = bvn * asr / 2.0 + norm_cdf(h1) * norm_cdf(hk)
    return max(0.0, min(1.0, bvn))


def _bs2002_call(S, K, t, r, b, sigma):
    """American call via the Bjerksund-Stensland 2002 two-region model."""
    # No early exercise if carry >= rate: American call == European call.
    if b >= r:
        return bsm_price(S, K, t, r, sigma, OptionType.CALL, b=b)

    v2 = sigma * sigma
    beta = (0.5 - b / v2) + math.sqrt((b / v2 - 0.5) ** 2 + 2.0 * r / v2)
    b_inf = beta / (beta - 1.0) * K
    b0 = max(K, r / (r - b) * K)

    t1 = 0.5 * (math.sqrt(5.0) - 1.0) * t

    h1 = -(b * t1 + 2.0 * sigma * math.sqrt(t1)) * K * K / ((b_inf - b0) * b0)
    h2 = -(b * t + 2.0 * sigma * math.sqrt(t)) * K * K / ((b_inf - b0) * b0)
    I1 = b0 + (b_inf - b0) * (1.0 - math.exp(h1))
    I2 = b0 + (b_inf - b0) * (1.0 - math.exp(h2))

    alpha1 = (I1 - K) * I1 ** (-beta)
    alpha2 = (I2 - K) * I2 ** (-beta)

    if S >= I2:
        return S - K  # immediate exercise region

    return (alpha2 * S ** beta
            - alpha2 * _phi(S, t1, beta, I2, I2, r, b, sigma)
            + _phi(S, t1, 1.0, I2, I2, r, b, sigma)
            - _phi(S, t1, 1.0, I1, I2, r, b, sigma)
            - K * _phi(S, t1, 0.0, I2, I2, r, b, sigma)
            + K * _phi(S, t1, 0.0, I1, I2, r, b, sigma)
            + alpha1 * _phi(S, t1, beta, I1, I2, r, b, sigma)
            - alpha1 * _psi(S, t, beta, I1, I2, I1, r, b, sigma, t1)
            + _psi(S, t, 1.0, I1, I2, I1, r, b, sigma, t1)
            - _psi(S, t, 1.0, K, I2, I1, r, b, sigma, t1)
            - K * _psi(S, t, 0.0, I1, I2, I1, r, b, sigma, t1)
            + K * _psi(S, t, 0.0, K, I2, I1, r, b, sigma, t1))
    # NOTE: _phi(S, t, gamma, h, I) and _psi(S, T2, gamma, h, I2, I1, t1).


def bjerksund_stensland(S, K, t, r, sigma, option_type=OptionType.CALL, b=None):
    """American option price via Bjerksund-Stensland (2002), closed form.

    Args mirror the rest of the engine. ``b`` is the cost of carry (defaults to
    ``r``); dividend yield q enters as b = r - q. American puts are priced via
    the exact put-call transformation P(S,K,r,b) = C(K,S,r-b,-b).
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        disc = math.exp(-r * t)
        payoff = max(fwd - K, 0.0) if ot is OptionType.CALL else max(K - fwd, 0.0)
        # Immediate-exercise intrinsic may exceed the discounted forward payoff.
        intrinsic = max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)
        return max(disc * payoff, intrinsic)

    if ot is OptionType.CALL:
        return _bs2002_call(S, K, t, r, b, sigma)
    # Put via transformation.
    return _bs2002_call(K, S, t, r - b, -b, sigma)


def bjerksund_stensland_greeks(S, K, t, r, sigma, option_type=OptionType.CALL,
                               b=None):
    """Greeks of the Bjerksund-Stensland American price by finite differences.

    The 2002 price is a closed form but its Greeks have no simple expression
    (the exercise boundary and the bivariate-normal term move with the inputs),
    so we central-difference the price. Returns a dict with delta, gamma, vega,
    theta (per year, calendar), and rho.

    Bumps are chosen small relative to each input; because the BS2002 price is a
    smooth function of its arguments (away from t=0) central differences are
    accurate to a few basis points, plenty for hedging.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, r_=r, sigma_=sigma, b_=b):
        return bjerksund_stensland(S_, K, t_, r_, sigma_, ot, b=b_)

    base = px()

    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)

    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)

    # Calendar theta = -dPrice/d(t_expiry); guard against t - h <= 0.
    ht = min(1e-4, 0.5 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)

    # Rho: bump the rate and the carry together (the stock case b moves with r).
    hr = 1e-5
    rho = (px(r_=r + hr, b_=b + hr) - px(r_=r - hr, b_=b - hr)) / (2.0 * hr)

    return {"price": base, "delta": delta, "gamma": gamma,
            "vega": vega, "theta": theta, "rho": rho}


def early_exercise_premium(S, K, t, r, sigma, option_type=OptionType.CALL,
                           b=None):
    """Decompose the American price into European value + early-exercise premium.

    Returns a dict with ``american`` (Bjerksund-Stensland), ``european`` (BSM),
    and ``premium`` = american - european, the extra value from the right to
    exercise early. The premium is non-negative and is (near) zero for an
    American call with no dividends (``b >= r``), where early exercise is never
    optimal.
    """
    from .bsm import price as bsm_price
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    american = bjerksund_stensland(S, K, t, r, sigma, ot, b=b)
    european = bsm_price(S, K, t, r, sigma, ot, b=b)
    premium = american - european
    return {"american": american, "european": european,
            "premium": max(premium, 0.0)}


def _bs1993_call(S, K, t, r, b, sigma):
    """American call via the Bjerksund-Stensland 1993 single-boundary model."""
    if b >= r:
        return bsm_price(S, K, t, r, sigma, OptionType.CALL, b=b)

    v2 = sigma * sigma
    beta = (0.5 - b / v2) + math.sqrt((b / v2 - 0.5) ** 2 + 2.0 * r / v2)
    b_inf = beta / (beta - 1.0) * K
    b0 = max(K, r / (r - b) * K)
    h = -(b * t + 2.0 * sigma * math.sqrt(t)) * (b0 / (b_inf - b0))
    trigger = b0 + (b_inf - b0) * (1.0 - math.exp(h))

    if S >= trigger:
        return S - K
    alpha = (trigger - K) * trigger ** (-beta)
    return (alpha * S ** beta
            - alpha * _phi(S, t, beta, trigger, trigger, r, b, sigma)
            + _phi(S, t, 1.0, trigger, trigger, r, b, sigma)
            - _phi(S, t, 1.0, K, trigger, r, b, sigma)
            - K * _phi(S, t, 0.0, trigger, trigger, r, b, sigma)
            + K * _phi(S, t, 0.0, K, trigger, r, b, sigma))


def bjerksund_stensland_1993(S, K, t, r, sigma, option_type=OptionType.CALL, b=None):
    """American option price via Bjerksund-Stensland (1993), single flat boundary.

    A simpler and slightly less accurate predecessor to the 2002 two-region
    model (:func:`bjerksund_stensland`): it uses one flat exercise boundary. Calls
    are priced directly; puts via the exact transformation
    ``P(S,K,r,b) = C(K,S,r-b,-b)``.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        disc = math.exp(-r * t)
        fwd = S * math.exp(b * t)
        payoff = max(fwd - K, 0.0) if ot is OptionType.CALL else max(K - fwd, 0.0)
        intrinsic = max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)
        return max(disc * payoff, intrinsic)
    if ot is OptionType.CALL:
        return _bs1993_call(S, K, t, r, b, sigma)
    return _bs1993_call(K, S, t, r - b, -b, sigma)
