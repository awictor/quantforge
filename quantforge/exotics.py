"""Closed-form prices for common exotic options.

Everything here is analytic (no simulation): binary/digital options,
single-barrier options via the Reiner-Rubinstein formulas, and geometric-
average Asian options. All use the generalized cost-of-carry ``b`` (defaults
to ``r``) so they cover stocks, dividends, futures and FX like the vanilla
engine.

Cross-checks enforced by the test suite:
  * cash-or-nothing + K * asset-parity relationships,
  * in-out parity: knock-in + knock-out = vanilla,
  * geometric Asian reduces to a BSM price with adjusted vol and carry.
"""

import math
from enum import Enum

from .mathfns import norm_cdf
from .bsm import OptionType, _coerce_type, _validate, price as bsm_price


# --------------------------------------------------------------------------
# Digital / binary options
# --------------------------------------------------------------------------
def cash_or_nothing(S, K, t, r, sigma, option_type=OptionType.CALL, b=None, cash=1.0):
    """Pays ``cash`` if the option finishes in the money, else 0.

    Call pays when S_T > K; put pays when S_T < K.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        itm = fwd > K if ot is OptionType.CALL else fwd < K
        return disc * cash if itm else 0.0
    vsqrt = sigma * math.sqrt(t)
    d2 = (math.log(S / K) + (b - 0.5 * sigma * sigma) * t) / vsqrt
    if ot is OptionType.CALL:
        return cash * disc * norm_cdf(d2)
    return cash * disc * norm_cdf(-d2)


def asset_or_nothing(S, K, t, r, sigma, option_type=OptionType.CALL, b=None):
    """Pays the asset value S_T if in the money, else 0."""
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        itm = fwd > K if ot is OptionType.CALL else fwd < K
        return S * math.exp((b - r) * t) if itm else 0.0
    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log(S / K) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    carry = math.exp((b - r) * t)
    if ot is OptionType.CALL:
        return S * carry * norm_cdf(d1)
    return S * carry * norm_cdf(-d1)


# --------------------------------------------------------------------------
# Single-barrier options (Reiner-Rubinstein / Merton)
# --------------------------------------------------------------------------
class Barrier(str, Enum):
    DOWN_IN = "down-in"
    DOWN_OUT = "down-out"
    UP_IN = "up-in"
    UP_OUT = "up-out"


def barrier_option(S, K, H, t, r, sigma, option_type=OptionType.CALL,
                   barrier=Barrier.DOWN_OUT, b=None, rebate=0.0):
    """Price a single-barrier option with an optional cash rebate.

    Args:
        H: barrier level.
        barrier: one of the four Barrier kinds.
        rebate: cash paid if the option is knocked out (out types) or never
            knocked in (in types), paid at expiry.

    Implements the standard Reiner-Rubinstein decomposition. Validated in the
    suite against in-out parity (knock-in + knock-out = vanilla + rebate term).
    """
    ot = _coerce_type(option_type)
    barrier = Barrier(barrier)
    _validate(S, K, t, sigma)
    if H <= 0:
        raise ValueError("barrier H must be positive")
    if b is None:
        b = r

    mu = (b - 0.5 * sigma * sigma) / (sigma * sigma)
    lam = math.sqrt(mu * mu + 2.0 * r / (sigma * sigma))
    vsqrt = sigma * math.sqrt(t)
    carry = math.exp((b - r) * t)
    disc = math.exp(-r * t)

    # phi = +1 for call, -1 for put ; eta = +1 for down, -1 for up
    phi = 1.0 if ot is OptionType.CALL else -1.0
    eta = 1.0 if barrier in (Barrier.DOWN_IN, Barrier.DOWN_OUT) else -1.0

    x1 = math.log(S / K) / vsqrt + (1.0 + mu) * vsqrt
    x2 = math.log(S / H) / vsqrt + (1.0 + mu) * vsqrt
    y1 = math.log(H * H / (S * K)) / vsqrt + (1.0 + mu) * vsqrt
    y2 = math.log(H / S) / vsqrt + (1.0 + mu) * vsqrt

    HS_2mu = (H / S) ** (2.0 * mu)
    HS_2mu2 = (H / S) ** (2.0 * mu + 2.0)

    A = (phi * S * carry * norm_cdf(phi * x1)
         - phi * K * disc * norm_cdf(phi * x1 - phi * vsqrt))
    B = (phi * S * carry * norm_cdf(phi * x2)
         - phi * K * disc * norm_cdf(phi * x2 - phi * vsqrt))
    C = (phi * S * carry * HS_2mu2 * norm_cdf(eta * y1)
         - phi * K * disc * HS_2mu * norm_cdf(eta * y1 - eta * vsqrt))
    D = (phi * S * carry * HS_2mu2 * norm_cdf(eta * y2)
         - phi * K * disc * HS_2mu * norm_cdf(eta * y2 - eta * vsqrt))

    # Rebate value: paid at expiry.
    z = math.log(H / S) / vsqrt + lam * vsqrt
    E = rebate * disc * (norm_cdf(eta * x2 - eta * vsqrt)
                         - HS_2mu * norm_cdf(eta * y2 - eta * vsqrt))
    F = rebate * ((H / S) ** (mu + lam) * norm_cdf(eta * z)
                  + (H / S) ** (mu - lam) * norm_cdf(eta * z - 2.0 * eta * lam * vsqrt))

    in_the_money_strike = K > H  # relative position of strike vs barrier

    # Assemble per Reiner-Rubinstein tables. Cases split on K vs H.
    call = ot is OptionType.CALL
    down = eta == 1.0

    if call:
        if down:  # down barriers, call
            if barrier in (Barrier.DOWN_IN,):
                val = (C + E) if in_the_money_strike else (A - B + D + E)
            else:  # DOWN_OUT
                val = (A - C + F) if in_the_money_strike else (B - D + F)
        else:  # up barriers, call
            if barrier in (Barrier.UP_IN,):
                val = (A + E) if in_the_money_strike else (B - C + D + E)
            else:  # UP_OUT
                val = F if in_the_money_strike else (A - B + C - D + F)
    else:  # put
        if down:  # down barriers, put
            if barrier in (Barrier.DOWN_IN,):
                val = (B - C + D + E) if in_the_money_strike else (A + E)
            else:  # DOWN_OUT
                val = (A - B + C - D + F) if in_the_money_strike else F
        else:  # up barriers, put
            if barrier in (Barrier.UP_IN,):
                val = (A - B + D + E) if in_the_money_strike else (C + E)
            else:  # UP_OUT
                val = (B - D + F) if in_the_money_strike else (A - C + F)
    return val


# --------------------------------------------------------------------------
# Geometric-average Asian option (closed form)
# --------------------------------------------------------------------------
def geometric_asian(S, K, t, r, sigma, option_type=OptionType.CALL, b=None):
    """Continuously-monitored geometric-average-price Asian option.

    The geometric average of a lognormal is itself lognormal, so the price is a
    BSM price with adjusted volatility and carry:

        sigma_A = sigma / sqrt(3)
        b_A     = 0.5 * (b - sigma^2 / 6)

    (Kemna-Vorst). This gives an exact closed form and is a standard control
    variate for the arithmetic-average Asian priced by simulation.
    """
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    sigma_a = sigma / math.sqrt(3.0)
    b_a = 0.5 * (b - sigma * sigma / 6.0)
    return bsm_price(S, K, t, r, sigma_a, option_type, b=b_a)


# --------------------------------------------------------------------------
# Arithmetic-average Asian option (Turnbull-Wakeman moment matching)
# --------------------------------------------------------------------------
def arithmetic_asian(S, K, t, r, sigma, option_type=OptionType.CALL, b=None):
    """Continuously-monitored arithmetic-average-price Asian (Turnbull-Wakeman).

    The arithmetic average of a lognormal is not lognormal, so there is no exact
    closed form. Turnbull-Wakeman (1991) matches the first two moments of the
    average to a lognormal and prices with a Black-Scholes-style formula on the
    average's forward. Fast and accurate for typical vols; agrees with the
    arithmetic-Asian Monte Carlo (:func:`quantforge.arithmetic_asian_mc`) to a
    few basis points. Averaging runs over the full life ``[0, t]``.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    # First moment of the continuous arithmetic average over [0, t].
    if abs(b) > 1e-12:
        M1 = S * (math.exp(b * t) - 1.0) / (b * t)
    else:
        M1 = S  # b -> 0 limit

    if t == 0 or sigma == 0:
        disc = math.exp(-r * t)
        payoff = max(M1 - K, 0.0) if ot is OptionType.CALL else max(K - M1, 0.0)
        return disc * payoff

    # Second moment (Turnbull-Wakeman).
    v2 = sigma * sigma
    if abs(b) > 1e-12:
        term1 = (2.0 * S * S * math.exp((2.0 * b + v2) * t)) / ((b + v2) * (2.0 * b + v2) * t * t)
        term2 = (2.0 * S * S / (b * t * t)) * (1.0 / (2.0 * b + v2)
                                               - math.exp(b * t) / (b + v2))
        M2 = term1 + term2
    else:
        M2 = (2.0 * S * S / (v2 * v2 * t * t)) * (math.exp(v2 * t) - 1.0 - v2 * t)

    # Match to a lognormal: effective vol and carry from the two moments.
    sigma_a = math.sqrt(math.log(M2 / (M1 * M1)) / t)
    b_a = math.log(M1 / S) / t
    return bsm_price(S, K, t, r, sigma_a, ot, b=b_a)


# --------------------------------------------------------------------------
# One-touch / no-touch binaries (touch options)
# --------------------------------------------------------------------------
def one_touch(S, H, t, r, sigma, b=None, cash=1.0, payoff_at_hit=True):
    """One-touch binary: pays ``cash`` if the spot ever reaches barrier ``H``.

    A continuously-monitored American digital. ``payoff_at_hit=True`` pays the
    cash immediately when the barrier is touched (the FX-market convention);
    ``False`` defers the payment to expiry. Works for an upper barrier
    (``H > S``) or a lower barrier (``H < S``); the direction is inferred.

    Uses the standard Rubinstein-Reiner touch formulas.
    """
    _validate(S, H, t, sigma)
    if H <= 0:
        raise ValueError("barrier H must be positive")
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        return 0.0  # cannot touch in zero time unless already there

    # eta = -1 for an up-barrier (H > S), +1 for a down-barrier (H < S).
    eta = 1.0 if H < S else -1.0
    vsqrt = sigma * math.sqrt(t)
    mu = (b - 0.5 * sigma * sigma) / (sigma * sigma)

    if payoff_at_hit:
        # Pays cash immediately on touch: lambda uses the full discount rate.
        lam = math.sqrt(mu * mu + 2.0 * r / (sigma * sigma))
        z = math.log(H / S) / vsqrt + lam * vsqrt
        HS = H / S
        return cash * (HS ** (mu + lam) * norm_cdf(eta * z)
                       + HS ** (mu - lam) * norm_cdf(eta * z - 2.0 * eta * lam * vsqrt))
    else:
        # Pays cash at expiry if touched: discount the barrier hit probability.
        # P(hit) = N(eta(-a + m)/v) + (H/S)^{2 mu} N(eta(-a - m)/v), with
        # a = ln(H/S), m = (b - sigma^2/2) t, eta = +1 up / -1 down barrier.
        disc = math.exp(-r * t)
        a = math.log(H / S)
        m = (b - 0.5 * sigma * sigma) * t
        # eta_dir = +1 for up-barrier (H > S), -1 for down-barrier.
        eta_dir = 1.0 if H > S else -1.0
        HS = H / S
        p_hit = (norm_cdf(eta_dir * (-a + m) / vsqrt)
                 + HS ** (2.0 * mu) * norm_cdf(eta_dir * (-a - m) / vsqrt))
        return cash * disc * p_hit


def no_touch(S, H, t, r, sigma, b=None, cash=1.0):
    """No-touch binary: pays ``cash`` at expiry if the barrier is never reached.

    Complementary to :func:`one_touch` with payment at expiry:
    ``no_touch = cash * e^{-rt} - one_touch(payoff_at_hit=False)``.
    """
    _validate(S, H, t, sigma)
    if b is None:
        b = r
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        return cash * disc  # never touches in zero time
    hit = one_touch(S, H, t, r, sigma, b=b, cash=cash, payoff_at_hit=False)
    return cash * disc - hit


# --------------------------------------------------------------------------
# Barrier option Greeks (finite differences on the closed form)
# --------------------------------------------------------------------------
def barrier_greeks(S, K, H, t, r, sigma, option_type=OptionType.CALL,
                   barrier=Barrier.DOWN_OUT, b=None, rebate=0.0):
    """Greeks of a single-barrier option by central finite differences.

    The Reiner-Rubinstein price is a closed form, but its Greeks are messy and
    change character across the barrier, so we central-difference the price.
    Returns a dict with delta, gamma, vega, and theta (calendar, per year).

    Near the barrier the true delta/gamma are large and discontinuous; the
    finite differences there are indicative rather than exact - use a spot bump
    well away from ``H`` when a smooth number is needed.
    """
    ot = _coerce_type(option_type)
    barrier = Barrier(barrier)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        return barrier_option(S_, K, H, t_, r, sigma_, ot, barrier, b=b,
                              rebate=rebate)

    base = px()
    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)

    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)

    ht = min(1e-4, 0.5 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)

    return {"price": base, "delta": delta, "gamma": gamma,
            "vega": vega, "theta": theta}


# --------------------------------------------------------------------------
# Asian option Greeks (finite differences on the closed forms)
# --------------------------------------------------------------------------
def asian_greeks(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                 average="geometric"):
    """Greeks of an Asian option by central finite differences.

    ``average`` selects the closed form to differentiate: "geometric"
    (Kemna-Vorst, exact) or "arithmetic" (Turnbull-Wakeman moment match).
    Returns a dict with delta, gamma, vega, and theta (calendar, per year).
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if average == "geometric":
        pricer = geometric_asian
    elif average == "arithmetic":
        pricer = arithmetic_asian
    else:
        raise ValueError("average must be 'geometric' or 'arithmetic'")

    def px(S_=S, t_=t, sigma_=sigma):
        return pricer(S_, K, t_, r, sigma_, ot, b=b)

    base = px()
    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)

    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)

    ht = min(1e-4, 0.5 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)

    return {"price": base, "delta": delta, "gamma": gamma,
            "vega": vega, "theta": theta}
