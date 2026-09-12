"""Displaced-diffusion (shifted lognormal) option pricing.

Rubinstein's displaced-diffusion model assumes ``S + shift`` is lognormal
rather than ``S`` itself. The shift lets the underlying reach values as low as
``-shift`` (so negative rates/prices are allowed) and produces a monotone skew:
a positive shift flattens the smile toward normal-like (Bachelier) behavior,
while ``shift = 0`` recovers Black-Scholes.

Pricing is exact: it is a Black-Scholes price on the displaced variables
``S' = S + shift`` and ``K' = K + shift`` with a rescaled volatility
``sigma' = sigma * S / (S + shift)`` chosen so the at-the-money instantaneous
volatility matches ``sigma`` (a common calibration convention).
"""

import math

from .bsm import price as bsm_price, OptionType, _coerce_type


def displaced_diffusion_price(S, K, t, r, sigma, shift=0.0,
                              option_type=OptionType.CALL, b=None) -> float:
    """Price a European option under the displaced-diffusion model.

    Args:
        shift: the displacement added to spot and strike. ``shift = 0`` is
            Black-Scholes; larger positive shifts push toward normal-model
            behavior and allow the underlying to fall below zero (down to
            ``-shift``).
        b: cost of carry (defaults to r).

    The payoff is unchanged (``max(S_T - K, 0)`` etc.); only the diffusion is
    displaced, so the price equals a BSM price on ``S + shift`` / ``K + shift``.
    """
    ot = _coerce_type(option_type)
    if t < 0:
        raise ValueError("t must be non-negative")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if b is None:
        b = r
    # Displaced diffusion permits negative spot/strike (down to -shift); only
    # the shifted variables must be positive.
    if S + shift <= 0 or K + shift <= 0:
        raise ValueError("S + shift and K + shift must be positive")

    S_shift = S + shift
    K_shift = K + shift
    # Rescale vol so the ATM instantaneous vol of the displaced process matches
    # sigma at the current spot: sigma_displaced = sigma * S / (S + shift).
    sigma_d = sigma * S / S_shift

    # The displaced forward must equal the true forward plus the (carried)
    # shift, so use a carry that reproduces F' = S*e^{bt} + shift on S_shift.
    fwd = S * math.exp(b * t)
    fwd_shift = fwd + shift
    b_d = math.log(fwd_shift / S_shift) / t if t > 0 else b

    return bsm_price(S_shift, K_shift, t, r, sigma_d, ot, b=b_d)


def displaced_diffusion_greeks(S, K, t, r, sigma, shift=0.0,
                               option_type=OptionType.CALL, b=None):
    """Greeks of a displaced-diffusion option by central finite differences.

    Differentiates :func:`displaced_diffusion_price` for ``delta`` (dV/dS),
    ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), and ``theta`` (calendar decay).
    At ``shift = 0`` these reduce to the vanilla Black-Scholes Greeks; a positive
    shift flattens the smile toward normal-model behaviour. Returns a dict with
    ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    if t < 0:
        raise ValueError("t must be non-negative")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if b is None:
        b = r
    if S + shift <= 0 or K + shift <= 0:
        raise ValueError("S + shift and K + shift must be positive")

    def px(S_=S, t_=t, sigma_=sigma):
        return displaced_diffusion_price(S_, K, t_, r, sigma_, shift, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t) if t > 0 else 1e-4
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def displaced_diffusion_smile(S, strikes, t, r, sigma, shift=0.0, b=None):
    """The Black-Scholes implied-vol smile a displaced-diffusion model produces.

    Prices a European call at each strike under the displaced diffusion, then
    inverts each price to its Black-Scholes implied volatility, returning
    ``(log_moneyness, vol)`` pairs sorted by strike (log-moneyness on the forward
    ``F = S e^{b t}``). Because ``sigma`` is calibrated to the ATM instantaneous
    vol, the smile passes near ``sigma`` at the money; a positive ``shift`` makes
    the process partly normal, producing a downward skew (steeper for larger
    shift), while ``shift = 0`` returns a flat Black-Scholes smile.
    """
    from .implied import implied_volatility
    if b is None:
        b = r
    F = S * math.exp(b * t)
    out = []
    for K in sorted(strikes):
        c = displaced_diffusion_price(S, K, t, r, sigma, shift,
                                      OptionType.CALL, b=b)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=b)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out


def displaced_implied_shift(S, t, r, quotes, b=None,
                            shift_lo=None, shift_hi=None):
    """Calibrate the displacement that reproduces an observed vol skew.

    Displaced diffusion has one skew knob, the ``shift``: a positive shift
    lowers the low-strike wing relative to the high-strike wing (a downward
    skew), while ``shift = 0`` is flat Black-Scholes. Given a set of Black-Scholes
    implied-vol quotes ``quotes = [(K, iv), ...]`` this finds the single shift
    whose displaced-diffusion smile best fits them.

    The at-the-money volatility is not a free skew knob here, but it is not
    fixed blindly either: for every trial shift the model's local ``sigma`` is
    re-solved so the displaced smile reproduces the ATM quote (the one whose
    strike is closest to the forward ``F = S e^{b t}``) exactly. That decouples
    level from skew, so the shift is driven purely by the off-ATM quotes and the
    fit is not biased by the local-vs-implied vol convention.

    Returns ``(shift, sigma_atm, rmse)`` where ``sigma_atm`` is the local vol at
    the fitted shift and ``rmse`` is the root-mean-square implied-vol error
    across the quotes. Minimises the squared vol error over the shift by
    golden-section search on ``[shift_lo, shift_hi]`` (defaults scale with spot:
    ``[-0.9 S, 20 S]``, staying above the ``-shift`` floor).
    """
    from .implied import implied_volatility

    if b is None:
        b = r
    quotes = [(float(K), float(iv)) for K, iv in quotes]
    if len(quotes) < 2:
        raise ValueError("need at least two (K, iv) quotes to fit a skew")
    if any(iv <= 0 for _, iv in quotes):
        raise ValueError("implied vols must be positive")

    F = S * math.exp(b * t)
    # The ATM anchor: quote nearest the forward, matched exactly at every shift.
    K_atm, iv_atm = min(quotes, key=lambda q: abs(q[0] - F))

    if shift_lo is None:
        shift_lo = -0.9 * S
    if shift_hi is None:
        shift_hi = 20.0 * S

    def sigma_for_shift(shift):
        # Solve the local sigma whose displaced price at K_atm inverts to iv_atm.
        # The map sigma -> displaced ATM implied vol is monotone, so bisect.
        lo, hi = 1e-6, 5.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            px = displaced_diffusion_price(S, K_atm, t, r, mid, shift,
                                           OptionType.CALL, b=b)
            iv = implied_volatility(px, S, K_atm, t, r, OptionType.CALL, b=b)
            if iv < iv_atm:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def sse(shift):
        if S + shift <= 0 or K_atm + shift <= 0:
            return float("inf"), 0.0
        sigma_loc = sigma_for_shift(shift)
        total = 0.0
        for K, iv_mkt in quotes:
            if K + shift <= 0:
                return float("inf"), sigma_loc
            px = displaced_diffusion_price(S, K, t, r, sigma_loc, shift,
                                           OptionType.CALL, b=b)
            try:
                iv_mod = implied_volatility(px, S, K, t, r,
                                            OptionType.CALL, b=b)
            except ValueError:
                return float("inf"), sigma_loc
            total += (iv_mod - iv_mkt) ** 2
        return total, sigma_loc

    # Golden-section search for the minimiser on [shift_lo, shift_hi].
    invphi = (math.sqrt(5.0) - 1.0) / 2.0
    a, bb = shift_lo, shift_hi
    c = bb - invphi * (bb - a)
    d = a + invphi * (bb - a)
    fc, fd = sse(c)[0], sse(d)[0]
    for _ in range(200):
        if bb - a < 1e-8 * max(1.0, S):
            break
        if fc < fd:
            bb, d, fd = d, c, fc
            c = bb - invphi * (bb - a)
            fc = sse(c)[0]
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (bb - a)
            fd = sse(d)[0]
    shift = 0.5 * (a + bb)
    total, sigma_atm = sse(shift)
    rmse = math.sqrt(total / len(quotes))
    return shift, sigma_atm, rmse


def displaced_diffusion_implied_vol(target_price, S, K, t, r, shift=0.0,
                                    option_type=OptionType.CALL, b=None,
                                    tol=1e-10, max_iter=200):
    """Implied displaced-diffusion volatility from a market price.

    Inverts :func:`displaced_diffusion_price` for the ``sigma`` reproducing
    ``target_price`` by bisection (the price is monotone increasing in ``sigma``).
    With ``shift = 0`` this coincides with the Black-Scholes implied vol. Raises if
    the quote lies outside the attainable ``[intrinsic, forward]`` band.
    """
    ot = _coerce_type(option_type)
    if t <= 0:
        raise ValueError("cannot imply vol at or past expiry")

    def price(sig):
        return displaced_diffusion_price(S, K, t, r, sig, shift, ot, b=b)

    lo, hi = 1e-9, 10.0
    p_lo, p_hi = price(lo), price(hi)
    if not (p_lo - 1e-12 <= target_price <= p_hi + 1e-12):
        raise ValueError("target price outside the attainable range")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        pm = price(mid)
        if abs(pm - target_price) < tol:
            return mid
        if pm < target_price:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
