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

from .mathfns import norm_cdf, norm_pdf
from .bsm import (
    OptionType, _coerce_type, _validate, price as bsm_price,
    delta as bsm_delta, gamma as bsm_gamma, vega as bsm_vega,
    theta as bsm_theta,
)


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


def range_binary(S, K_low, K_high, t, r, sigma, b=None, cash=1.0):
    """Range binary (double digital): pays ``cash`` iff ``K_low <= S_T <= K_high``.

    A bet that the terminal price lands inside a corridor at expiry (no path
    monitoring). It is exactly the difference of two cash-or-nothing calls
    struck at the two levels: paying when ``S_T > K_low`` but not when
    ``S_T > K_high``. With ``d2(K) = (log(S/K) + (b - sigma^2/2) t)/(sigma sqrt t)``,

        price = cash e^{-rt} (N(d2(K_low)) - N(d2(K_high))).
    """
    _validate(S, K_low, t, sigma)
    if K_high <= K_low:
        raise ValueError("K_high must exceed K_low")
    if b is None:
        b = r
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        return disc * cash if K_low <= fwd <= K_high else 0.0
    vsqrt = sigma * math.sqrt(t)
    d2_lo = (math.log(S / K_low) + (b - 0.5 * sigma * sigma) * t) / vsqrt
    d2_hi = (math.log(S / K_high) + (b - 0.5 * sigma * sigma) * t) / vsqrt
    return cash * disc * (norm_cdf(d2_lo) - norm_cdf(d2_hi))


def supershare(S, K_low, K_high, t, r, sigma, b=None):
    """Supershare option: pays ``S_T / K_low`` iff ``K_low <= S_T <= K_high``.

    Introduced by Hakansson (1976) as a building block for mutual-fund payoffs.
    It is a scaled difference of two asset-or-nothing calls struck at the two
    levels, so with ``d1(K) = (log(S/K) + (b + sigma^2/2) t)/(sigma sqrt t)`` and
    carry ``e^{(b-r)t}``,

        price = (S/K_low) e^{(b-r)t} (N(d1(K_low)) - N(d1(K_high))).
    """
    _validate(S, K_low, t, sigma)
    if K_high <= K_low:
        raise ValueError("K_high must exceed K_low")
    if b is None:
        b = r
    carry = math.exp((b - r) * t)
    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        return (S / K_low) * carry if K_low <= fwd <= K_high else 0.0
    vsqrt = sigma * math.sqrt(t)
    d1_lo = (math.log(S / K_low) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    d1_hi = (math.log(S / K_high) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    return (S / K_low) * carry * (norm_cdf(d1_lo) - norm_cdf(d1_hi))


def range_binary_greeks(S, K_low, K_high, t, r, sigma, b=None, cash=1.0):
    """Greeks of a range binary (double digital).

    ``delta`` and ``gamma`` are analytic. With ``d2(K)`` the digital exponent
    and ``e^{-rt}`` the discount, the corridor value is
    ``cash e^{-rt} (N(d2_lo) - N(d2_hi))``, so differentiating in spot,

        delta = cash e^{-rt} (phi(d2_lo) - phi(d2_hi)) / (S sigma sqrt(t))
        gamma = -cash e^{-rt} / (S^2 sigma sqrt(t))
                * ((phi(d2_lo) d2_lo - phi(d2_hi) d2_hi) / (sigma sqrt(t))
                   + phi(d2_lo) - phi(d2_hi)).

    ``vega`` and ``theta`` (calendar decay) are central finite differences of
    :func:`range_binary`. Returns a dict with ``price`` and those fields. The
    delta changes sign across the middle of the corridor and gamma is large near
    either edge as expiry approaches (double-sided pin risk).
    """
    _validate(S, K_low, t, sigma)
    if K_high <= K_low:
        raise ValueError("K_high must exceed K_low")
    if b is None:
        b = r
    disc = math.exp(-r * t)
    vsqrt = sigma * math.sqrt(t)
    d2_lo = (math.log(S / K_low) + (b - 0.5 * sigma * sigma) * t) / vsqrt
    d2_hi = (math.log(S / K_high) + (b - 0.5 * sigma * sigma) * t) / vsqrt
    plo, phi = norm_pdf(d2_lo), norm_pdf(d2_hi)
    price = cash * disc * (norm_cdf(d2_lo) - norm_cdf(d2_hi))
    delta = cash * disc * (plo - phi) / (S * vsqrt)
    # d(phi(d2)/(S vsqrt))/dS with d(d2)/dS = 1/(S vsqrt):
    #   d/dS [phi(d2)] = -d2 phi(d2)/(S vsqrt); plus the 1/S from the prefactor.
    gamma = -cash * disc / (S * S * vsqrt) * (
        (plo * d2_lo - phi * d2_hi) / vsqrt + plo - phi)

    def px(sigma_=sigma, t_=t):
        return range_binary(S, K_low, K_high, t_, r, sigma_, b=b, cash=cash)

    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": price, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def supershare_greeks(S, K_low, K_high, t, r, sigma, b=None):
    """Greeks of a supershare option by central finite differences of
    :func:`supershare`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2), ``vega``
    (dV/dsigma), ``theta`` (calendar decay). Returns a dict with ``price`` and
    those fields.
    """
    _validate(S, K_low, t, sigma)
    if K_high <= K_low:
        raise ValueError("K_high must exceed K_low")
    if b is None:
        b = r

    def px(S_=S, sigma_=sigma, t_=t):
        return supershare(S_, K_low, K_high, t_, r, sigma_, b=b)

    base = px()
    carry = math.exp((b - r) * t)
    vsqrt = sigma * math.sqrt(t)
    d1_lo = (math.log(S / K_low) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    d1_hi = (math.log(S / K_high) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    plo, phi = norm_pdf(d1_lo), norm_pdf(d1_hi)
    dN = norm_cdf(d1_lo) - norm_cdf(d1_hi)
    dphi = plo - phi
    A = carry / K_low
    # price = A S (N(d1_lo) - N(d1_hi)); d(d1)/dS = 1/(S vsqrt).
    delta = A * (dN + dphi / vsqrt)
    # gamma = d(delta)/dS, using d(phi(d1))/dS = -d1 phi(d1)/(S vsqrt).
    gamma = A / (S * vsqrt) * (
        dphi - (d1_lo * plo - d1_hi * phi) / vsqrt)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


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


def geometric_asian_greeks(S, K, t, r, sigma, option_type=OptionType.CALL,
                           b=None):
    """Greeks of a continuously-monitored geometric-average Asian option.

    The Kemna-Vorst price is exactly a Black-Scholes price with the adjusted
    volatility ``sigma_A = sigma / sqrt(3)`` and carry ``b_A = (b - sigma^2/6)/2``,
    so the spot ``S`` enters only through the BSM price at ``(sigma_A, b_A)``:
    ``delta`` and ``gamma`` are the exact BSM Greeks evaluated there (no finite
    difference). ``vega``, ``theta``, and ``rho`` do depend on ``sigma``/``t``/
    ``r`` through the adjusted parameters, so they are taken as central finite
    differences of the exact closed form. Returns a dict with ``price``,
    ``delta``, ``gamma``, ``vega``, ``theta``, ``rho``.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    sigma_a = sigma / math.sqrt(3.0)
    b_a = 0.5 * (b - sigma * sigma / 6.0)
    price = bsm_price(S, K, t, r, sigma_a, ot, b=b_a)
    # S enters only via the BSM price at the adjusted params -> exact Greeks.
    delta = bsm_delta(S, K, t, r, sigma_a, ot, b=b_a)
    gamma = bsm_gamma(S, K, t, r, sigma_a, b=b_a)

    def px(ss=sigma, tt=t, rr=r):
        return geometric_asian(S, K, tt, rr, ss, ot, b=b)

    hv = 1e-4
    vega = (px(ss=sigma + hv) - px(ss=sigma - hv)) / (2.0 * hv)
    ht = 1e-4
    # Theta is -dV/dt (value decays as time passes).
    theta = -(px(tt=t + ht) - px(tt=t - ht)) / (2.0 * ht)
    hr = 1e-5
    rho = (px(rr=r + hr) - px(rr=r - hr)) / (2.0 * hr)
    return {"price": price, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta, "rho": rho}


def _geom_fixing_times(t, n):
    """Equally-spaced fixing times ``t*i/n`` for ``i = 1..n`` (last at expiry)."""
    if n < 1:
        raise ValueError("n_fixings must be >= 1")
    return [t * i / n for i in range(1, n + 1)]


def average_strike_arithmetic_asian(S, t, r, sigma, n_fixings=None,
                                    fixing_times=None,
                                    option_type=OptionType.CALL, b=None):
    """Average-strike (floating-strike) discrete arithmetic Asian option.

    The strike is the realized arithmetic average: a call pays
    ``max(S_T - A, 0)`` and a put ``max(A - S_T, 0)``, where
    ``A = (1/n) sum_i S_{t_i}``. This is an exchange option between the terminal
    price ``S_T`` (exactly lognormal) and the average ``A`` (matched to a
    lognormal by its first two moments, Levy 1992). The cross-moment
    ``E[S_T A]`` is exact,

        E[S_T A] = (S^2/n) sum_i exp(b (t + t_i) + sigma^2 min(t, t_i))
                 = (S^2/n) sum_i exp(b (t + t_i) + sigma^2 t_i),

    so the log-space covariance is ``rho_log = log(E[S_T A]/(E[S_T] E[A]))`` and
    the spread variance is ``Var = sigma^2 t + V_A - 2 rho_log`` with
    ``V_A = log(M2/M1^2)``. A call is then
    ``e^{-rt} (E[S_T] N(d1) - E[A] N(d2))`` in the usual Margrabe form.

    Provide either ``n_fixings`` (equally-spaced dates ``t*i/n``, last at
    expiry) or an explicit ``fixing_times`` sequence in ``(0, t]``.
    """
    ot = _coerce_type(option_type)
    _validate(S, 1.0, t, sigma)
    if b is None:
        b = r
    if fixing_times is None:
        if n_fixings is None:
            raise ValueError("provide n_fixings or fixing_times")
        times = _geom_fixing_times(t, n_fixings)
    else:
        times = [float(x) for x in fixing_times]
        if not times:
            raise ValueError("fixing_times must be non-empty")
        if any(x <= 0.0 or x > t + 1e-12 for x in times):
            raise ValueError("fixing_times must lie in (0, t]")
    n = len(times)
    disc = math.exp(-r * t)
    v2 = sigma * sigma
    # Arithmetic-average moments (Levy).
    M1 = (S / n) * sum(math.exp(b * ti) for ti in times)
    dbl = 0.0
    for ti in times:
        for tj in times:
            dbl += math.exp(b * (ti + tj) + v2 * (ti if ti < tj else tj))
    M2 = (S * S / (n * n)) * dbl
    E_ST = S * math.exp(b * t)
    E_A = M1
    # Exact cross-moment E[S_T A]; here min(t, t_i) = t_i since t_i <= t.
    E_ST_A = (S * S / n) * sum(math.exp(b * (t + ti) + v2 * ti) for ti in times)
    V_A = math.log(M2 / (M1 * M1))
    rho_log = math.log(E_ST_A / (E_ST * E_A))
    var = v2 * t + V_A - 2.0 * rho_log
    if var <= 0.0:
        payoff = (max(E_ST - E_A, 0.0) if ot is OptionType.CALL
                  else max(E_A - E_ST, 0.0))
        return disc * payoff
    sd = math.sqrt(var)
    d1 = (math.log(E_ST / E_A) + 0.5 * var) / sd
    d2 = d1 - sd
    if ot is OptionType.CALL:
        return disc * (E_ST * norm_cdf(d1) - E_A * norm_cdf(d2))
    return disc * (E_A * norm_cdf(-d2) - E_ST * norm_cdf(-d1))


def average_strike_arithmetic_asian_greeks(S, t, r, sigma, n_fixings=None,
                                           fixing_times=None,
                                           option_type=OptionType.CALL, b=None):
    """Greeks of an average-strike arithmetic Asian by central finite
    differences of :func:`average_strike_arithmetic_asian`: ``delta`` (dV/dS),
    ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay).
    Returns a dict with ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, 1.0, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        ft = (None if fixing_times is None
              else [x * t_ / t for x in fixing_times])
        return average_strike_arithmetic_asian(S_, t_, r, sigma_, n_fixings,
                                               ft, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def average_strike_geometric_asian(S, t, r, sigma, n_fixings=None,
                                   fixing_times=None,
                                   option_type=OptionType.CALL, b=None):
    """Average-strike (floating-strike) discrete geometric Asian option (exact).

    The strike is the realized geometric average rather than a fixed level:
    a call pays ``max(S_T - G, 0)`` and a put ``max(G - S_T, 0)``, where
    ``G = (prod_i S_{t_i})^{1/n}`` is the geometric average over the fixing
    dates. The terminal price ``S_T`` and the average ``G`` are jointly
    lognormal, so this is an exchange option between two lognormal assets and
    has an exact Margrabe-style closed form.

    With ``E[S_T] = S e^{bt}``, ``E[G] = exp(m_G + v_G/2)`` (the lognormal
    average forward), and the variance of ``log S_T - log G``

        Var = sigma^2 t + v_G - 2 sigma^2 mean(t_i),
        v_G = (sigma^2/n^2) sum_ij min(t_i, t_j),

    a call is ``e^{-rt} (E[S_T] N(d1) - E[G] N(d2))`` with
    ``d1 = (log(E[S_T]/E[G]) + Var/2)/sqrt(Var)`` and ``d2 = d1 - sqrt(Var)``.

    Provide either ``n_fixings`` (equally-spaced dates ``t*i/n``, last at
    expiry) or an explicit ``fixing_times`` sequence in ``(0, t]``.
    """
    ot = _coerce_type(option_type)
    _validate(S, 1.0, t, sigma)  # K unused; validate S, t, sigma
    if b is None:
        b = r
    if fixing_times is None:
        if n_fixings is None:
            raise ValueError("provide n_fixings or fixing_times")
        times = _geom_fixing_times(t, n_fixings)
    else:
        times = [float(x) for x in fixing_times]
        if not times:
            raise ValueError("fixing_times must be non-empty")
        if any(x <= 0.0 or x > t + 1e-12 for x in times):
            raise ValueError("fixing_times must lie in (0, t]")
    n = len(times)
    disc = math.exp(-r * t)
    v2 = sigma * sigma
    mean_t = sum(times) / n
    # Lognormal parameters of the geometric average G.
    m_G = math.log(S) + (b - 0.5 * v2) * mean_t
    dbl = 0.0
    for ti in times:
        for tj in times:
            dbl += ti if ti < tj else tj
    v_G = v2 * dbl / (n * n)
    E_ST = S * math.exp(b * t)
    E_G = math.exp(m_G + 0.5 * v_G)
    # cov(log S_T, log G) = (sigma^2/n) sum_i min(T, t_i) = sigma^2 mean(t_i).
    cov = v2 * mean_t
    var = v2 * t + v_G - 2.0 * cov
    if var <= 0.0:
        payoff = (max(E_ST - E_G, 0.0) if ot is OptionType.CALL
                  else max(E_G - E_ST, 0.0))
        return disc * payoff
    sd = math.sqrt(var)
    d1 = (math.log(E_ST / E_G) + 0.5 * var) / sd
    d2 = d1 - sd
    if ot is OptionType.CALL:
        return disc * (E_ST * norm_cdf(d1) - E_G * norm_cdf(d2))
    return disc * (E_G * norm_cdf(-d2) - E_ST * norm_cdf(-d1))


def average_strike_geometric_asian_greeks(S, t, r, sigma, n_fixings=None,
                                          fixing_times=None,
                                          option_type=OptionType.CALL, b=None):
    """Greeks of an average-strike geometric Asian by central finite
    differences of :func:`average_strike_geometric_asian`: ``delta`` (dV/dS),
    ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay).
    Returns a dict with ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, 1.0, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        # Scale the schedule with t so a calendar bump preserves the fixing
        # shape and stays within (0, t_].
        ft = (None if fixing_times is None
              else [x * t_ / t for x in fixing_times])
        return average_strike_geometric_asian(S_, t_, r, sigma_, n_fixings,
                                              ft, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def discrete_geometric_asian(S, K, t, r, sigma, n_fixings=None,
                             fixing_times=None, option_type=OptionType.CALL,
                             b=None):
    """Discretely-monitored geometric-average-price Asian option (exact).

    The geometric average ``G = (prod_i S_{t_i})^{1/n}`` over the monitoring
    dates ``t_i`` is lognormal, because ``log G`` is a linear combination of the
    jointly-Gaussian log-prices. With ``m = log S + (b - sigma^2/2) * mean(t_i)``
    and ``v = (sigma^2 / n^2) * sum_i sum_j min(t_i, t_j)``, ``log G`` is
    ``Normal(m, v)`` and the price is a Black-Scholes-style closed form on the
    forward ``F = exp(m + v/2)`` discounted at ``r``:

        d1 = (m + v - log K) / sqrt(v),  d2 = d1 - sqrt(v)
        call = e^{-rt} (F N(d1) - K N(d2)).

    Provide either ``n_fixings`` (equally-spaced dates ``t*i/n``, last at expiry)
    or an explicit ``fixing_times`` sequence in ``(0, t]``. A single fixing at
    ``t`` recovers the vanilla Black-Scholes price; as ``n_fixings -> infinity``
    the price converges to the continuous Kemna-Vorst :func:`geometric_asian`.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if fixing_times is None:
        if n_fixings is None:
            raise ValueError("provide n_fixings or fixing_times")
        times = _geom_fixing_times(t, n_fixings)
    else:
        times = [float(x) for x in fixing_times]
        if not times:
            raise ValueError("fixing_times must be non-empty")
        if any(x <= 0.0 or x > t + 1e-12 for x in times):
            raise ValueError("fixing_times must lie in (0, t]")
    n = len(times)
    disc = math.exp(-r * t)
    mean_t = sum(times) / n
    m = math.log(S) + (b - 0.5 * sigma * sigma) * mean_t
    # v = (sigma^2 / n^2) * sum_i sum_j min(t_i, t_j).
    dbl = 0.0
    for ti in times:
        for tj in times:
            dbl += ti if ti < tj else tj
    v = sigma * sigma * dbl / (n * n)
    if v <= 0.0:
        fwd = math.exp(m)
        payoff = max(fwd - K, 0.0) if ot is OptionType.CALL else max(K - fwd, 0.0)
        return disc * payoff
    fwd = math.exp(m + 0.5 * v)
    sd = math.sqrt(v)
    d1 = (m + v - math.log(K)) / sd
    d2 = d1 - sd
    if ot is OptionType.CALL:
        return disc * (fwd * norm_cdf(d1) - K * norm_cdf(d2))
    return disc * (K * norm_cdf(-d2) - fwd * norm_cdf(-d1))


def discrete_geometric_asian_greeks(S, K, t, r, sigma, n_fixings=None,
                                    fixing_times=None,
                                    option_type=OptionType.CALL, b=None):
    """Greeks of a discrete geometric-average Asian option by central finite
    differences of :func:`discrete_geometric_asian`: ``delta`` (dV/dS),
    ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay). When
    ``fixing_times`` is given it is held fixed; with ``n_fixings`` the equally-
    spaced grid rescales with ``t`` (matching the continuous convention).
    Returns a dict with ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        return discrete_geometric_asian(S_, K, t_, r, sigma_, n_fixings,
                                        fixing_times, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


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


def seasoned_geometric_asian(S, K, t, r, sigma, observed_prices, n_total,
                             remaining_fixing_times=None,
                             option_type=OptionType.CALL, b=None):
    """Seasoned (in-progress) discrete geometric-average Asian option (exact).

    Prices a geometric Asian partway through its averaging window, when some
    fixings have already been observed. With ``n_total`` fixings in all, of
    which ``observed_prices`` are already fixed and ``k`` remain at future
    times ``remaining_fixing_times`` (in ``(0, t]``, measured from now), the
    final average ``G = (prod_{i=1}^{n} S_{t_i})^{1/n}`` is still lognormal: the
    observed factors contribute a known constant ``A = sum log(S_obs)`` and the
    ``k`` future log-prices are jointly Gaussian. Hence ``log G`` is
    ``Normal(m, v)`` with

        m = (A + k log S + (b - sigma^2/2) sum_j tau_j) / n
        v = (sigma^2 / n^2) sum_i sum_j min(tau_i, tau_j),

    and the price is a Black-Scholes-style closed form on ``F = exp(m + v/2)``
    discounted at ``r``. When no fixings are observed this reduces exactly to
    :func:`discrete_geometric_asian`. When all ``n_total`` fixings are observed
    the average is known and the payoff is deterministic.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    obs = [float(x) for x in observed_prices]
    if any(x <= 0.0 for x in obs):
        raise ValueError("observed_prices must be positive")
    m_obs = len(obs)
    if n_total < 1:
        raise ValueError("n_total must be >= 1")
    if m_obs > n_total:
        raise ValueError("more observed fixings than n_total")
    k = n_total - m_obs
    if remaining_fixing_times is None:
        # Default: the remaining fixings are the last k of an equally-spaced
        # grid over [0, t] (so tau_j = t*(m_obs + j)/n_total for j=1..k, but
        # re-based to now they span (0, t]). Use an evenly-spaced future grid.
        times = [t * i / k for i in range(1, k + 1)] if k > 0 else []
    else:
        times = [float(x) for x in remaining_fixing_times]
        if len(times) != k:
            raise ValueError("remaining_fixing_times must have length n_total - len(observed_prices)")
        if any(x <= 0.0 or x > t + 1e-12 for x in times):
            raise ValueError("remaining_fixing_times must lie in (0, t]")
    disc = math.exp(-r * t)
    A = sum(math.log(x) for x in obs)
    n = n_total
    if k == 0:
        G = math.exp(A / n)
        payoff = max(G - K, 0.0) if ot is OptionType.CALL else max(K - G, 0.0)
        return disc * payoff
    m_mean = (A + k * math.log(S) + (b - 0.5 * sigma * sigma) * sum(times)) / n
    dbl = 0.0
    for ti in times:
        for tj in times:
            dbl += ti if ti < tj else tj
    v = sigma * sigma * dbl / (n * n)
    if v <= 0.0:
        fwd = math.exp(m_mean)
        payoff = max(fwd - K, 0.0) if ot is OptionType.CALL else max(K - fwd, 0.0)
        return disc * payoff
    fwd = math.exp(m_mean + 0.5 * v)
    sd = math.sqrt(v)
    d1 = (m_mean + v - math.log(K)) / sd
    d2 = d1 - sd
    if ot is OptionType.CALL:
        return disc * (fwd * norm_cdf(d1) - K * norm_cdf(d2))
    return disc * (K * norm_cdf(-d2) - fwd * norm_cdf(-d1))


def seasoned_geometric_asian_greeks(S, K, t, r, sigma, observed_prices,
                                    n_total, remaining_fixing_times=None,
                                    option_type=OptionType.CALL, b=None):
    """Greeks of a seasoned geometric Asian by central finite differences of
    :func:`seasoned_geometric_asian`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2),
    ``vega`` (dV/dsigma), ``theta`` (calendar decay). The observed fixings and
    the remaining schedule are held fixed. Returns a dict with ``price`` and
    those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        # Scale the remaining schedule with t so a calendar bump keeps the
        # fixing shape (and stays within (0, t_]); matches the n_fixings-
        # rescales-with-t convention of the fresh discrete greeks.
        rem = (None if remaining_fixing_times is None
               else [x * t_ / t for x in remaining_fixing_times])
        return seasoned_geometric_asian(S_, K, t_, r, sigma_, observed_prices,
                                        n_total, rem, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def seasoned_arithmetic_asian(S, K, t, r, sigma, observed_prices, n_total,
                              remaining_fixing_times=None,
                              option_type=OptionType.CALL, b=None):
    """Seasoned (in-progress) discrete arithmetic-average Asian (Levy match).

    Prices an arithmetic Asian partway through its averaging window, when some
    fixings are already observed. The final average is
    ``A = (Q + sum_j S_{tau_j}) / n`` where ``Q = sum(observed_prices)`` is a
    known constant and the ``k`` remaining prices are lognormal. A call payoff
    ``max(A - K, 0) = (1/n) max(sum_j S_{tau_j} - (n K - Q), 0)`` is therefore an
    arithmetic-average option on the *remaining* fixings with the shifted strike
    ``K' = n K - Q``, scaled by ``1/n``. The remaining sum's first two moments
    are exact,

        m1 = S sum_j exp(b tau_j)
        m2 = S^2 sum_ij exp(b (tau_i + tau_j) + sigma^2 min(tau_i, tau_j)),

    and Levy (1992) matches them to a lognormal priced by Black-Scholes.

    Special cases handled exactly: if ``K' <= 0`` the call is always in the
    money and worth ``e^{-rt} (E[A] - K)`` (the put is worthless), and vice
    versa. With no observations this reduces to :func:`discrete_arithmetic_asian`;
    with all fixings observed the payoff is the deterministic arithmetic
    intrinsic.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    obs = [float(x) for x in observed_prices]
    if any(x <= 0.0 for x in obs):
        raise ValueError("observed_prices must be positive")
    m_obs = len(obs)
    if n_total < 1:
        raise ValueError("n_total must be >= 1")
    if m_obs > n_total:
        raise ValueError("more observed fixings than n_total")
    k = n_total - m_obs
    if remaining_fixing_times is None:
        times = [t * i / k for i in range(1, k + 1)] if k > 0 else []
    else:
        times = [float(x) for x in remaining_fixing_times]
        if len(times) != k:
            raise ValueError("remaining_fixing_times must have length n_total - len(observed_prices)")
        if any(x <= 0.0 or x > t + 1e-12 for x in times):
            raise ValueError("remaining_fixing_times must lie in (0, t]")
    disc = math.exp(-r * t)
    n = n_total
    Q = sum(obs)
    if k == 0:
        A = Q / n
        payoff = max(A - K, 0.0) if ot is OptionType.CALL else max(K - A, 0.0)
        return disc * payoff
    v2 = sigma * sigma
    m1 = S * sum(math.exp(b * tj) for tj in times)
    dbl = 0.0
    for ti in times:
        for tj in times:
            dbl += math.exp(b * (ti + tj) + v2 * (ti if ti < tj else tj))
    m2 = S * S * dbl
    Kp = n * K - Q  # shifted strike on the remaining sum
    # Expected average for the always-in-the-money branches.
    EA = (Q + m1) / n
    if Kp <= 0.0:
        # Remaining sum is always >= Kp (=0 or negative): call always ITM.
        if ot is OptionType.CALL:
            return disc * (EA - K)
        return 0.0
    V = math.log(m2 / (m1 * m1))
    if V <= 0.0:
        payoff = max(m1 - Kp, 0.0) if ot is OptionType.CALL else max(Kp - m1, 0.0)
        return disc * payoff / n
    sd = math.sqrt(V)
    d1 = (math.log(m1 / Kp) + 0.5 * V) / sd
    d2 = d1 - sd
    if ot is OptionType.CALL:
        return disc * (m1 * norm_cdf(d1) - Kp * norm_cdf(d2)) / n
    return disc * (Kp * norm_cdf(-d2) - m1 * norm_cdf(-d1)) / n


def seasoned_arithmetic_asian_greeks(S, K, t, r, sigma, observed_prices,
                                     n_total, remaining_fixing_times=None,
                                     option_type=OptionType.CALL, b=None):
    """Greeks of a seasoned arithmetic Asian by central finite differences of
    :func:`seasoned_arithmetic_asian`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2),
    ``vega`` (dV/dsigma), ``theta`` (calendar decay). The observed fixings and
    the remaining schedule are held fixed. Returns a dict with ``price`` and
    those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        rem = (None if remaining_fixing_times is None
               else [x * t_ / t for x in remaining_fixing_times])
        return seasoned_arithmetic_asian(S_, K, t_, r, sigma_, observed_prices,
                                         n_total, rem, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def discrete_arithmetic_asian(S, K, t, r, sigma, n_fixings=None,
                              fixing_times=None, option_type=OptionType.CALL,
                              b=None):
    """Discretely-monitored arithmetic-average-price Asian option (Levy
    moment-matching approximation).

    The arithmetic average ``A = (1/n) sum_i S_{t_i}`` is not lognormal, but its
    first two moments over the fixing dates have exact closed forms:

        M1 = (S/n) sum_i exp(b t_i)
        M2 = (S^2/n^2) sum_i sum_j exp(b (t_i + t_j) + sigma^2 min(t_i, t_j)).

    Levy (1992) matches these to a lognormal and prices with a Black-Scholes
    formula on the average's forward ``M1`` and effective variance
    ``V = log(M2/M1^2)``:

        d1 = (log(M1/K) + V/2) / sqrt(V),  d2 = d1 - sqrt(V)
        call = e^{-rt} (M1 N(d1) - K N(d2)).

    Provide either ``n_fixings`` (equally-spaced dates ``t*i/n``, last at expiry)
    or an explicit ``fixing_times`` sequence in ``(0, t]``. A single fixing at
    ``t`` recovers the vanilla Black-Scholes price. The geometric-average
    Asian is an exact lower bound; this arithmetic price sits above it.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if fixing_times is None:
        if n_fixings is None:
            raise ValueError("provide n_fixings or fixing_times")
        times = _geom_fixing_times(t, n_fixings)
    else:
        times = [float(x) for x in fixing_times]
        if not times:
            raise ValueError("fixing_times must be non-empty")
        if any(x <= 0.0 or x > t + 1e-12 for x in times):
            raise ValueError("fixing_times must lie in (0, t]")
    n = len(times)
    disc = math.exp(-r * t)
    v2 = sigma * sigma
    M1 = (S / n) * sum(math.exp(b * ti) for ti in times)
    dbl = 0.0
    for ti in times:
        for tj in times:
            dbl += math.exp(b * (ti + tj) + v2 * (ti if ti < tj else tj))
    M2 = (S * S / (n * n)) * dbl
    V = math.log(M2 / (M1 * M1))
    if V <= 0.0:
        payoff = max(M1 - K, 0.0) if ot is OptionType.CALL else max(K - M1, 0.0)
        return disc * payoff
    sd = math.sqrt(V)
    d1 = (math.log(M1 / K) + 0.5 * V) / sd
    d2 = d1 - sd
    if ot is OptionType.CALL:
        return disc * (M1 * norm_cdf(d1) - K * norm_cdf(d2))
    return disc * (K * norm_cdf(-d2) - M1 * norm_cdf(-d1))


def discrete_arithmetic_asian_greeks(S, K, t, r, sigma, n_fixings=None,
                                     fixing_times=None,
                                     option_type=OptionType.CALL, b=None):
    """Greeks of a discrete arithmetic-average Asian option by central finite
    differences of :func:`discrete_arithmetic_asian`: ``delta`` (dV/dS),
    ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay).
    Returns a dict with ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        return discrete_arithmetic_asian(S_, K, t_, r, sigma_, n_fixings,
                                         fixing_times, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


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


def double_no_touch(S, L, U, t, r, sigma, b=None, cash=1.0, n_terms=200):
    """Double-no-touch: pays ``cash`` at expiry if spot stays inside ``(L, U)``.

    Continuously monitored: the option survives only if the spot never touches
    either the lower barrier ``L`` or the upper barrier ``U`` before expiry. The
    survival probability of driftful Brownian motion in a strip has the classic
    Fourier (eigenfunction) expansion; with ``x = ln(S/L)``, ``Z = ln(U/L)``,
    ``m = b - sigma^2/2`` and ``beta = m/sigma^2``,

        P(survive) = (2/Z) e^{-beta x - m^2 t / (2 sigma^2)}
            * sum_{n>=1} sin(k_n x) e^{-k_n^2 sigma^2 t / 2}
                         * k_n (1 - (-1)^n e^{beta Z}) / (beta^2 + k_n^2),

    with ``k_n = n pi / Z``. The value is ``cash e^{-rt} P(survive)``. As
    ``U -> infinity`` it collapses to the single lower :func:`no_touch`, and as
    ``L -> 0`` to the upper one. Requires ``L < S < U``.
    """
    _validate(S, L, t, sigma)
    if not (0.0 < L < S < U):
        raise ValueError("need 0 < L < S < U")
    if b is None:
        b = r
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        return cash * disc  # cannot touch in zero time while strictly inside

    x0 = math.log(S / L)
    Z = math.log(U / L)
    m = b - 0.5 * sigma * sigma
    beta = m / (sigma * sigma)
    var = sigma * sigma * t
    pref = (2.0 / Z) * math.exp(-beta * x0 - m * m * t / (2.0 * sigma * sigma))
    total = 0.0
    for n in range(1, n_terms + 1):
        kn = n * math.pi / Z
        total += (math.sin(kn * x0) * math.exp(-0.5 * kn * kn * var)
                  * kn * (1.0 - ((-1) ** n) * math.exp(beta * Z))
                  / (beta * beta + kn * kn))
    surv = pref * total
    # Clamp tiny numerical negatives/overshoots from the truncated series.
    surv = min(1.0, max(0.0, surv))
    return cash * disc * surv


def double_one_touch(S, L, U, t, r, sigma, b=None, cash=1.0, n_terms=200):
    """Double-one-touch: pays ``cash`` at expiry if spot touches ``L`` or ``U``.

    The expiry-settled complement of :func:`double_no_touch`:
    ``double_one_touch = cash e^{-rt} - double_no_touch``. Requires ``L < S < U``.
    """
    disc = math.exp(-r * t)
    dnt = double_no_touch(S, L, U, t, r, sigma, b=b, cash=cash, n_terms=n_terms)
    return cash * disc - dnt


def double_knock_out_call(S, K, L, U, t, r, sigma, b=None, delta1=0.0,
                          delta2=0.0, n_terms=10):
    """Ikeda-Kunitomo (1992) double-barrier knock-out call: payoff max(S_T - K, 0).

    Pays the vanilla call payoff only if the continuously-monitored spot stays
    strictly inside the (possibly exponentially curved) corridor bounded below by
    ``L e^{delta1 s}`` and above by ``U e^{delta2 s}`` over ``[0, t]``. With
    ``delta1 = delta2 = 0`` the barriers are flat at ``L`` and ``U``. The price is
    the Ikeda-Kunitomo image series

        C = S e^{(b-r)t} sum_n [ (U^n/L^n)^{mu1} (L^n/S)^{mu2} (N(d1)-N(d2))
                                 - (L^{n+1}/(U^n S))^{mu3} (N(d3)-N(d4)) ]
            - K e^{-rt} sum_n [ ... same with mu-2 and d-sigma sqrt(t) ... ],

    truncated at ``|n| <= n_terms`` (the series converges geometrically). Requires
    ``L < S < U`` and ``K < U`` for a non-trivial payoff.
    """
    _validate(S, K, t, sigma)
    if not (0.0 < L < S < U):
        raise ValueError("need 0 < L < S < U")
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        return max(S - K, 0.0) if L < S < U else 0.0

    vt = sigma * math.sqrt(t)
    V = sigma * sigma
    F = U  # upper barrier level in the payoff cap
    E = max(K, L)  # the call is worthless below max(K, L)
    disc_q = math.exp((b - r) * t)
    disc_r = math.exp(-r * t)

    sum1 = 0.0
    sum2 = 0.0
    for n in range(-n_terms, n_terms + 1):
        mu1 = 2.0 * (b - delta2 - n * (delta1 - delta2)) / V + 1.0
        mu2 = 2.0 * n * (delta1 - delta2) / V
        mu3 = 2.0 * (b - delta2 + n * (delta1 - delta2)) / V + 1.0
        d1 = (math.log(S * U ** (2 * n) / (E * L ** (2 * n)))
              + (b + 0.5 * V) * t) / vt
        d2 = (math.log(S * U ** (2 * n) / (F * L ** (2 * n)))
              + (b + 0.5 * V) * t) / vt
        d3 = (math.log(L ** (2 * n + 2) / (E * S * U ** (2 * n)))
              + (b + 0.5 * V) * t) / vt
        d4 = (math.log(L ** (2 * n + 2) / (F * S * U ** (2 * n)))
              + (b + 0.5 * V) * t) / vt
        term1 = ((U ** n / L ** n) ** mu1 * (L ** n / S) ** mu2
                 * (norm_cdf(d1) - norm_cdf(d2))
                 - (L ** (n + 1) / (U ** n * S)) ** mu3
                 * (norm_cdf(d3) - norm_cdf(d4)))
        term2 = ((U ** n / L ** n) ** (mu1 - 2.0) * (L ** n / S) ** mu2
                 * (norm_cdf(d1 - vt) - norm_cdf(d2 - vt))
                 - (L ** (n + 1) / (U ** n * S)) ** (mu3 - 2.0)
                 * (norm_cdf(d3 - vt) - norm_cdf(d4 - vt)))
        sum1 += term1
        sum2 += term2
    return max(S * disc_q * sum1 - K * disc_r * sum2, 0.0)


def double_knock_in_call(S, K, L, U, t, r, sigma, b=None, delta1=0.0,
                         delta2=0.0, n_terms=10):
    """Double-barrier knock-in call: pays the call only if a barrier is touched.

    The in-out complement of :func:`double_knock_out_call`: a knock-in and a
    knock-out with the same strike and corridor partition every path, so at expiry

        double_knock_in_call + double_knock_out_call = vanilla call.

    Priced as ``vanilla - double_knock_out_call`` with the vanilla evaluated on the
    same carry ``b``. Requires ``L < S < U``.
    """
    if b is None:
        b = r
    vanilla = bsm_price(S, K, t, r, sigma, OptionType.CALL, b=b)
    dko = double_knock_out_call(S, K, L, U, t, r, sigma, b=b, delta1=delta1,
                                delta2=delta2, n_terms=n_terms)
    return max(vanilla - dko, 0.0)


def double_knock_out_call_greeks(S, K, L, U, t, r, sigma, b=None, delta1=0.0,
                                 delta2=0.0, n_terms=10):
    """Greeks of an Ikeda-Kunitomo double knock-out call by finite differences.

    Central differences of :func:`double_knock_out_call` for ``delta`` (dV/dS),
    ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay
    ``-dV/dt``), and the two barrier sensitivities ``dV/dL`` and ``dV/dU``. The
    value is a knock-out, so more volatility raises the knock probability and
    ``vega < 0`` near the middle of the corridor, and widening either barrier
    raises the value (``dV/dL < 0``, ``dV/dU > 0``). Returns a dict with those
    fields.
    """
    _validate(S, K, t, sigma)
    if not (0.0 < L < S < U):
        raise ValueError("need 0 < L < S < U")
    if b is None:
        b = r

    def px(S_=S, sigma_=sigma, t_=t, L_=L, U_=U):
        return double_knock_out_call(S_, K, L_, U_, t_, r, sigma_, b=b,
                                     delta1=delta1, delta2=delta2,
                                     n_terms=n_terms)

    base = px()
    # Keep spot bumps well inside the corridor so S +/- h stays in (L, U).
    hS = min(1e-4 * S, 0.25 * (S - L), 0.25 * (U - S))
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    hL = min(1e-4 * S, 0.25 * (S - L))
    dV_dL = (px(L_=L + hL) - px(L_=L - hL)) / (2.0 * hL)
    hU = min(1e-4 * S, 0.25 * (U - S))
    dV_dU = (px(U_=U + hU) - px(U_=U - hU)) / (2.0 * hU)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta, "dV_dL": dV_dL, "dV_dU": dV_dU}


def double_knock_in_call_greeks(S, K, L, U, t, r, sigma, b=None, delta1=0.0,
                                delta2=0.0, n_terms=10):
    """Greeks of a double-barrier knock-in call by in-out parity.

    Differentiating ``knock_in = vanilla - knock_out`` term by term, the spot,
    vol and time Greeks are the vanilla Black-Scholes Greek minus the double
    knock-out Greek (:func:`double_knock_out_call_greeks`). The vanilla has no
    barrier dependence, so the knock-in's barrier sensitivities are the negatives
    of the knock-out's (widening the corridor lowers the knock-in). Returns the
    same dict layout: ``price``, ``delta``, ``gamma``, ``vega``, ``theta``,
    ``dV_dL``, ``dV_dU``.
    """
    if b is None:
        b = r
    ot = OptionType.CALL
    ko = double_knock_out_call_greeks(S, K, L, U, t, r, sigma, b=b,
                                      delta1=delta1, delta2=delta2,
                                      n_terms=n_terms)
    # Vanilla Greeks (theta by the same -dV/dt convention as the KO greeks).
    ht = min(1e-4, 0.25 * t)
    van_theta = -(bsm_price(S, K, t + ht, r, sigma, ot, b=b)
                  - bsm_price(S, K, t - ht, r, sigma, ot, b=b)) / (2.0 * ht)
    return {
        "price": bsm_price(S, K, t, r, sigma, ot, b=b) - ko["price"],
        "delta": bsm_delta(S, K, t, r, sigma, ot, b=b) - ko["delta"],
        "gamma": bsm_gamma(S, K, t, r, sigma, b=b) - ko["gamma"],
        "vega": bsm_vega(S, K, t, r, sigma, b=b) - ko["vega"],
        "theta": van_theta - ko["theta"],
        "dV_dL": -ko["dV_dL"],
        "dV_dU": -ko["dV_dU"],
    }


def double_no_touch_greeks(S, L, U, t, r, sigma, b=None, cash=1.0, n_terms=200):
    """Greeks of a double-no-touch by finite differences on the closed form.

    Central differences of :func:`double_no_touch` for ``delta`` (dV/dS),
    ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay
    ``-dV/dt``), and the two barrier sensitivities ``dV/dL`` and ``dV/dU``. A DNT
    is a bet on low realized range, so ``vega < 0`` (more vol -> more likely to
    knock) and widening either barrier raises the value (``dV/dL < 0`` since a
    lower ``L`` widens the band, ``dV/dU > 0``). Returns a dict with those fields.
    """
    _validate(S, L, t, sigma)
    if not (0.0 < L < S < U):
        raise ValueError("need 0 < L < S < U")
    if b is None:
        b = r

    def px(S_=S, sigma_=sigma, t_=t, L_=L, U_=U):
        return double_no_touch(S_, L_, U_, t_, r, sigma_, b=b, cash=cash,
                               n_terms=n_terms)

    base = px()
    # Keep spot bumps well inside the band so S +/- h stays in (L, U).
    hS = min(1e-4 * S, 0.25 * (S - L), 0.25 * (U - S))
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    hL = min(1e-4 * S, 0.25 * (S - L))
    dV_dL = (px(L_=L + hL) - px(L_=L - hL)) / (2.0 * hL)
    hU = min(1e-4 * S, 0.25 * (U - S))
    dV_dU = (px(U_=U + hU) - px(U_=U - hU)) / (2.0 * hU)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta, "dV_dL": dV_dL, "dV_dU": dV_dU}


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


# --------------------------------------------------------------------------
# Gap options (Reiner-Rubinstein): separate trigger and payoff strikes
# --------------------------------------------------------------------------
def gap_option(S, K_trigger, K_payoff, t, r, sigma, option_type=OptionType.CALL,
               b=None):
    """Gap option: pays off against ``K_payoff`` but is triggered by ``K_trigger``.

    A gap call pays ``S_T - K_payoff`` (which may be negative) whenever
    ``S_T > K_trigger``; a gap put pays ``K_payoff - S_T`` whenever
    ``S_T < K_trigger``. Setting the two strikes equal recovers the vanilla
    option. Closed form (Reiner-Rubinstein).
    """
    ot = _coerce_type(option_type)
    _validate(S, K_trigger, t, sigma)
    if K_payoff <= 0:
        raise ValueError("payoff strike must be positive")
    if b is None:
        b = r
    carry = math.exp((b - r) * t)
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        if ot is OptionType.CALL:
            return disc * ((fwd - K_payoff) if fwd > K_trigger else 0.0)
        return disc * ((K_payoff - fwd) if fwd < K_trigger else 0.0)
    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log(S / K_trigger) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    d2 = d1 - vsqrt
    if ot is OptionType.CALL:
        return S * carry * norm_cdf(d1) - K_payoff * disc * norm_cdf(d2)
    return K_payoff * disc * norm_cdf(-d2) - S * carry * norm_cdf(-d1)


def gap_option_greeks(S, K_trigger, K_payoff, t, r, sigma,
                      option_type=OptionType.CALL, b=None):
    """Greeks of a gap option (Reiner-Rubinstein).

    ``delta`` and ``gamma`` are analytic; ``vega`` and ``theta`` (calendar
    decay) are central finite differences. A gap option is an asset-or-nothing
    minus ``K_payoff`` cash-or-nothings, both triggered at ``K_trigger``, so
    with ``d1, d2`` at the trigger and the identity
    ``S carry phi(d1) = K_trigger disc phi(d2)`` the spot sensitivities collapse
    to (call)

        delta = carry N(d1) + disc (K_trigger - K_payoff) phi(d2)/(S sigma sqrt t)
        gamma = carry phi(d1)/(S sigma sqrt t)
                - disc (K_trigger - K_payoff) phi(d2) (d2/(sigma sqrt t) + 1)/(S^2 sigma sqrt t).

    Setting ``K_trigger = K_payoff`` recovers the vanilla Black-Scholes delta
    and gamma. Returns a dict with ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, K_trigger, t, sigma)
    if K_payoff <= 0:
        raise ValueError("payoff strike must be positive")
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        return gap_option(S_, K_trigger, K_payoff, t_, r, sigma_, ot, b=b)

    base = px()
    carry = math.exp((b - r) * t)
    disc = math.exp(-r * t)
    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log(S / K_trigger) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    d2 = d1 - vsqrt
    gap = K_trigger - K_payoff
    pdf2 = norm_pdf(d2)
    if ot is OptionType.CALL:
        delta = carry * norm_cdf(d1) + disc * gap * pdf2 / (S * vsqrt)
    else:
        delta = -carry * norm_cdf(-d1) + disc * gap * pdf2 / (S * vsqrt)
    # gamma is the same for calls and puts (the N(+/-d1) term contributes
    # carry phi(d1)/(S vsqrt) either way; the cash term's second derivative is
    # symmetric in the sign convention).
    gamma = (carry * norm_pdf(d1) / (S * vsqrt)
             - disc * gap * pdf2 * (d2 / vsqrt + 1.0) / (S * S * vsqrt))
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


# --------------------------------------------------------------------------
# Power options: payoff on S^power
# --------------------------------------------------------------------------
def power_option(S, K, t, r, sigma, power, option_type=OptionType.CALL, b=None):
    """Power option with payoff ``max(S_T^power - K, 0)`` (call) / ``max(K - S_T^power, 0)``.

    S_T^power is lognormal, so this has a closed form: an adjusted-drift,
    adjusted-vol Black-Scholes on the transformed underlying. ``power = 1``
    recovers the vanilla option.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if power <= 0:
        raise ValueError("power must be positive")
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        fwd_pow = (S * math.exp(b * t)) ** power
        payoff = max(fwd_pow - K, 0.0) if ot is OptionType.CALL else max(K - fwd_pow, 0.0)
        return disc * payoff
    # E[S_T^power] under the risk-neutral measure with carry b.
    vsqrt = sigma * math.sqrt(t)
    # Mean and vol of log(S_T^power) = power * log(S_T).
    mu = power * (math.log(S) + (b - 0.5 * sigma * sigma) * t)
    vol_p = power * vsqrt
    # Forward of S^power and its lognormal parameters.
    fwd = math.exp(mu + 0.5 * vol_p * vol_p)
    d1 = (mu + vol_p * vol_p - math.log(K)) / vol_p
    d2 = d1 - vol_p
    if ot is OptionType.CALL:
        return disc * (fwd * norm_cdf(d1) - K * norm_cdf(d2))
    return disc * (K * norm_cdf(-d2) - fwd * norm_cdf(-d1))


def power_option_greeks(S, K, t, r, sigma, power, option_type=OptionType.CALL,
                        b=None):
    """Greeks of a power option (payoff ``max(S_T^power - K, 0)``) by FD.

    Central finite differences of the closed-form :func:`power_option` for
    ``delta`` (dV/dS), ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), and ``theta``
    (calendar decay, ``-dV/dt``). At ``power = 1`` these reduce to the vanilla
    Black-Scholes Greeks. Returns a dict with ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if power <= 0:
        raise ValueError("power must be positive")

    def px(S_=S, t_=t, sigma_=sigma):
        return power_option(S_, K, t_, r, sigma_, power, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def powered_option(S, K, t, r, sigma, power, option_type=OptionType.CALL,
                   b=None):
    """Powered option: payoff ``max(S_T - K, 0)**power`` (call) or
    ``max(K - S_T, 0)**power`` (put), for a positive **integer** ``power``.

    Distinct from :func:`power_option` (whose payoff is ``max(S_T**power - K, 0)``):
    here the *option payoff itself* is raised to a power, so the payoff has a
    higher-order convexity in the terminal spot. Because the payoff is a
    polynomial in ``S_T`` on the exercise region, it decomposes by the binomial
    theorem into a sum of ``S_T**j`` truncated moments, each of which has a
    closed form (Esser 2003; Heynen-Kat 1996). ``power = 1`` recovers the
    vanilla Black-Scholes option.

    With ``F_j = E[S_T**j] = S**j exp(j b t + 0.5 j (j-1) sigma^2 t)`` and
    ``d_j = (ln(S/K) + (b + (j - 0.5) sigma^2) t) / (sigma sqrt(t))``, a call is
    ``disc * sum_j C(p,j) (-K)^{p-j} F_j N(d_j)`` and a put is
    ``disc * sum_j C(p,j) K^{p-j} (-1)^j F_j N(-d_j)``.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if not isinstance(power, int) or power < 1:
        raise ValueError("power must be a positive integer")
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        intrinsic = (max(fwd - K, 0.0) if ot is OptionType.CALL
                     else max(K - fwd, 0.0))
        return disc * intrinsic ** power
    vsqrt = sigma * math.sqrt(t)
    total = 0.0
    for j in range(power + 1):
        c = math.comb(power, j)
        fwd_j = S ** j * math.exp(j * b * t + 0.5 * j * (j - 1) * sigma * sigma * t)
        d_j = (math.log(S / K) + (b + (j - 0.5) * sigma * sigma) * t) / vsqrt
        if ot is OptionType.CALL:
            total += c * (-K) ** (power - j) * fwd_j * norm_cdf(d_j)
        else:
            total += c * K ** (power - j) * (-1) ** j * fwd_j * norm_cdf(-d_j)
    return disc * total


def powered_option_greeks(S, K, t, r, sigma, power,
                          option_type=OptionType.CALL, b=None):
    """Greeks of a powered option by central finite differences of
    :func:`powered_option`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2), ``vega``
    (dV/dsigma), ``theta`` (calendar decay, ``-dV/dt``). Returns a dict with
    ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if not isinstance(power, int) or power < 1:
        raise ValueError("power must be a positive integer")

    def px(S_=S, t_=t, sigma_=sigma):
        return powered_option(S_, K, t_, r, sigma_, power, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def barrier_rebate(S, H, t, r, sigma, knock="out", b=None, cash=1.0,
                   payoff_at_hit=True):
    """Standalone rebate cashflow attached to a barrier.

    A **knock-out rebate** pays ``cash`` if the barrier ``H`` is breached (the
    consolation for the option knocking out); a **knock-in rebate** pays ``cash``
    at expiry if the barrier is *never* breached (the option failed to knock in).

    ``payoff_at_hit`` (knock-out only) pays on touch vs at expiry. This reuses
    the touch-option machinery: a knock-out rebate is a one-touch, a knock-in
    rebate is a no-touch.
    """
    knock = str(knock).lower()
    if knock in ("out", "knock-out", "ko"):
        return one_touch(S, H, t, r, sigma, b=b, cash=cash,
                         payoff_at_hit=payoff_at_hit)
    if knock in ("in", "knock-in", "ki"):
        # A knock-in rebate pays only if the barrier is never hit -> no-touch.
        return no_touch(S, H, t, r, sigma, b=b, cash=cash)
    raise ValueError("knock must be 'out' or 'in'")


def digital_greeks(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                   cash=1.0):
    """Delta and gamma of a cash-or-nothing digital by finite differences.

    Returns a dict with price, delta, and gamma. Near the strike as expiry
    approaches, the digital's delta spikes (and gamma flips sign across the
    strike) -- the "pin risk" that makes digitals hard to hedge and motivates
    the call-spread over-hedge in :mod:`quantforge.overhedge`.
    """
    _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    def px(S_):
        return cash_or_nothing(S_, K, t, r, sigma, option_type, b=b, cash=cash)

    base = px(S)
    h = 1e-3 * S
    up, dn = px(S + h), px(S - h)
    delta = (up - dn) / (2.0 * h)
    gamma = (up - 2.0 * base + dn) / (h * h)
    return {"price": base, "delta": delta, "gamma": gamma}
