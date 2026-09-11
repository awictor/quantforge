"""Extendible options (Longstaff 1990), closed form.

A holder-extendible call gives its owner, at the first expiry ``t1``, three
choices: exercise against the initial strike ``K1`` for ``S_{t1} - K1``, let the
option lapse, or pay a fee ``A`` to extend the life to ``T2`` with a (usually
different) strike ``K2``. The holder takes whichever is worth most:

    payoff(t1) = max(S_{t1} - K1, C(S_{t1}, K2, T2 - t1) - A, 0).

Longstaff (1990) prices this in closed form. The value is a vanilla call to the
first expiry plus a correction for the extension region, expressed with
bivariate normals that couple ``t1`` to the extended expiry ``T2``.

The extension is exercised on an interval of terminal spots ``[I_low, I_high]``:
below ``I_low`` the option lapses, above ``I_high`` immediate exercise beats
extending. The two boundaries solve ``C(I, K2, T2 - t1) - A = 0`` (lower) and
``C(I, K2, T2 - t1) - A = I - K1`` (upper).
"""

import math

from .bsm import call_price, put_price, _validate, OptionType, _coerce_type
from .american import _bivariate_normal
from .mathfns import norm_cdf


def _extend_boundaries(K1, K2, tau, r, sigma, b, A):
    """Solve the two terminal-spot boundaries of the extension region.

    ``I_low``: extended call value equals the fee (``C(I, K2, tau) = A``).
    ``I_high``: extending is no better than exercising now
    (``C(I, K2, tau) - A = I - K1``).
    """
    def cval(x):
        return call_price(x, K2, tau, r, sigma, b=b)

    # Lower boundary: C(I, K2, tau) = A, C increasing in spot.
    lo, hi = 1e-8, max(K1, K2)
    while cval(hi) < A:
        hi *= 2.0
        if hi > 1e12:
            break
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cval(mid) > A:
            hi = mid
        else:
            lo = mid
    I_low = 0.5 * (lo + hi)

    # Upper boundary: C(I, K2, tau) - A = I - K1  <=>  g(I) = C - A - I + K1 = 0.
    # g is decreasing in I for large I (call delta < 1), positive near I_low.
    def g(x):
        return cval(x) - A - x + K1

    lo2, hi2 = I_low, max(2.0 * K1, 2.0 * K2, I_low * 2.0)
    while g(hi2) > 0.0:
        hi2 *= 2.0
        if hi2 > 1e12:
            break
    for _ in range(200):
        mid = 0.5 * (lo2 + hi2)
        if g(mid) > 0.0:
            lo2 = mid
        else:
            hi2 = mid
    I_high = 0.5 * (lo2 + hi2)
    return I_low, I_high


def holder_extendible_call(S, K1, K2, t1, T2, r, sigma, A, b=None) -> float:
    """Holder-extendible call (Longstaff 1990), closed form.

    Args:
        K1: strike applying at the first expiry ``t1``.
        K2: strike of the extended option (life ``T2``).
        t1: first expiry (years); T2: extended expiry (years, ``> t1``).
        A: fee paid at ``t1`` to extend.

    Reduces to a vanilla call struck at ``K1`` expiring at ``t1`` when extending
    is never worthwhile (very large ``A``).
    """
    ot = OptionType.CALL
    if b is None:
        b = r
    _validate(S, K1, t1, sigma)
    _validate(S, K2, T2, sigma)
    if T2 <= t1:
        raise ValueError("require T2 > t1")
    if A < 0.0:
        raise ValueError("fee A must be non-negative")

    tau = T2 - t1
    I_low, I_high = _extend_boundaries(K1, K2, tau, r, sigma, b, A)

    # If the fee is so large that extending never beats exercising or lapsing,
    # the extension strip collapses (I_low >= I_high) and the option is simply a
    # vanilla call struck at K1 expiring at t1.
    if I_low >= I_high:
        return call_price(S, K1, t1, r, sigma, b=b)

    # The terminal payoff at t1 splits into three spot regions:
    #   S_{t1} < I_low        -> 0 (lapse)
    #   I_low <= S <= I_high  -> C(S, K2, tau) - A (extend)
    #   S > I_high            -> S - K1 (exercise now)
    # Value each region as a discounted expectation.
    v1 = sigma * math.sqrt(t1)
    st2 = sigma * math.sqrt(T2)
    rho = math.sqrt(t1 / T2)
    carry_T2 = math.exp((b - r) * T2)
    disc_T2 = math.exp(-r * T2)
    carry_t1 = math.exp((b - r) * t1)
    disc_t1 = math.exp(-r * t1)
    mu = b + 0.5 * sigma * sigma

    def z(level):
        return (math.log(S / level) + mu * t1) / v1

    z_lo, z_hi = z(I_low), z(I_high)
    y_K2 = (math.log(S / K2) + mu * T2) / st2

    # Extension strip [I_low, I_high]: extended-call value received, minus fee A.
    ext_asset = S * carry_T2 * (_bivariate_normal(z_lo, y_K2, rho)
                                - _bivariate_normal(z_hi, y_K2, rho))
    ext_cash = K2 * disc_T2 * (_bivariate_normal(z_lo - v1, y_K2 - st2, rho)
                               - _bivariate_normal(z_hi - v1, y_K2 - st2, rho))
    strip_cash = disc_t1 * (norm_cdf(z_lo - v1) - norm_cdf(z_hi - v1))
    extension_region = ext_asset - ext_cash - A * strip_cash

    # Exercise region S_{t1} > I_high: discounted E[(S - K1) 1{S > I_high}].
    exercise_region = (S * carry_t1 * norm_cdf(z_hi)
                       - K1 * disc_t1 * norm_cdf(z_hi - v1))

    return extension_region + exercise_region


def writer_extendible_put(S, K1, K2, t1, T2, r, sigma, b=None) -> float:
    """Writer-extendible put (Longstaff 1990), closed form.

    At the first expiry ``t1`` the put is exercised if it finishes in the money
    (``S_{t1} < K1``, paying ``K1 - S_{t1}``); otherwise the writer's obligation
    is automatically extended to ``T2`` as a put struck at ``K2`` (no fee). The
    terminal-``t1`` payoff is therefore

        payoff(t1) = (K1 - S_{t1})           if S_{t1} < K1
                   = P(S_{t1}, K2, T2 - t1)  if S_{t1} >= K1.

    The value is a vanilla put to ``t1`` plus the extended-put value collected on
    ``S_{t1} >= K1``, expressed with bivariate normals coupling ``t1`` and
    ``T2`` (correlation ``rho = sqrt(t1/T2)``):

        W = p(S, K1, t1)
            + K2 e^{-r T2} M(z2, -y2; -rho) - S e^{(b-r)T2} M(z1, -y1; -rho),

    with ``z2, y2`` the ``d2``-type arguments at ``K1`` (over ``t1``) and ``K2``
    (over ``T2``) and ``z1 = z2 + sigma sqrt(t1)``, ``y1 = y2 + sigma sqrt(T2)``.
    """
    if b is None:
        b = r
    _validate(S, K1, t1, sigma)
    _validate(S, K2, T2, sigma)
    if T2 <= t1:
        raise ValueError("require T2 > t1")

    v1 = sigma * math.sqrt(t1)
    st2 = sigma * math.sqrt(T2)
    rho = math.sqrt(t1 / T2)
    mu = b - 0.5 * sigma * sigma
    z2 = (math.log(S / K1) + mu * t1) / v1
    y2 = (math.log(S / K2) + mu * T2) / st2
    z1 = z2 + v1
    y1 = y2 + st2
    base = put_price(S, K1, t1, r, sigma, b=b)
    extension = (K2 * math.exp(-r * T2) * _bivariate_normal(z2, -y2, -rho)
                 - S * math.exp((b - r) * T2) * _bivariate_normal(z1, -y1, -rho))
    return base + extension


def writer_extendible_put_greeks(S, K1, K2, t1, T2, r, sigma, b=None):
    """Greeks of a writer-extendible put by central finite differences of
    :func:`writer_extendible_put`: ``delta``, ``gamma``, ``vega``, ``theta``
    (calendar decay, both expiries shrinking together). Returns a dict with
    ``price`` and those fields.
    """
    if b is None:
        b = r

    def px(S_=S, sigma_=sigma, shift=0.0):
        return writer_extendible_put(S_, K1, K2, t1 - shift, T2 - shift, r,
                                     sigma_, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t1)
    theta = -(px(shift=ht) - px(shift=-ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def holder_extendible_call_greeks(S, K1, K2, t1, T2, r, sigma, A, b=None):
    """Greeks of a holder-extendible call by central finite differences of
    :func:`holder_extendible_call`: ``delta``, ``gamma``, ``vega``, ``theta``
    (calendar decay, both expiries shrinking together). Returns a dict with
    ``price`` and those fields.
    """
    if b is None:
        b = r

    def px(S_=S, sigma_=sigma, shift=0.0):
        return holder_extendible_call(S_, K1, K2, t1 - shift, T2 - shift, r,
                                      sigma_, A, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t1)
    theta = -(px(shift=ht) - px(shift=-ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}
