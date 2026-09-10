"""CBOE VIX-style fair volatility index from an option chain or smile.

The VIX methodology computes a model-free 30-day (or arbitrary-tenor) implied
variance from a strip of out-of-the-money options and reports its square root,
scaled by 100:

    sigma^2 = (2 e^{r t} / t) * sum_i (dK_i / K_i^2) Q(K_i)
              - (1 / t) * (F / K0 - 1)^2,
    VIX = 100 * sqrt(sigma^2),

where ``F`` is the implied forward, ``K0`` the first strike below ``F``, ``Q``
the OTM option price at each strike (put below ``K0``, call above, average at
``K0``), and ``dK_i`` the central strike spacing. This module implements the
exact discrete CBOE formula (:func:`vix_from_chain`) and a convenience that
builds the chain from a smile (:func:`vix_from_smile`). Pure standard library.
"""

import math


def vix_from_chain(strikes, q_prices, F, t, r):
    """CBOE variance/VIX from a discrete OTM option chain.

    Args:
        strikes: increasing list of strikes.
        q_prices: the OTM option mid-price at each strike (put below the forward,
            call above; at the money use the average of the two).
        F: implied forward. t: tenor in years. r: risk-free rate.

    Returns ``(variance, vix)`` -- the annualized fair variance and
    ``100 * sqrt(variance)``.
    """
    n = len(strikes)
    if n != len(q_prices) or n < 3:
        raise ValueError("need at least three matching (strike, price) points")
    # K0 = first strike at or below F.
    k0 = strikes[0]
    for K in strikes:
        if K <= F:
            k0 = K
        else:
            break

    growth = math.exp(r * t)
    total = 0.0
    for i in range(n):
        if i == 0:
            dK = strikes[1] - strikes[0]
        elif i == n - 1:
            dK = strikes[n - 1] - strikes[n - 2]
        else:
            dK = 0.5 * (strikes[i + 1] - strikes[i - 1])
        total += (dK / (strikes[i] * strikes[i])) * q_prices[i]

    var = (2.0 * growth / t) * total - (1.0 / t) * (F / k0 - 1.0) ** 2
    return var, 100.0 * math.sqrt(max(var, 0.0))


def vix_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=201, width=6.0):
    """VIX-style fair index built from a smile ``vol_fn(K)``.

    Samples an OTM chain around the forward, prices each option at its smile vol
    with Black-Scholes, and applies :func:`vix_from_chain`. A flat smile returns
    ``VIX ~= 100 * sigma``.
    """
    from .bsm import call_price, put_price

    if t <= 0:
        raise ValueError("t must be positive")
    F = S0 * math.exp((r - q) * t)
    atm_vol = vol_fn(F)
    sd = atm_vol * math.sqrt(t)
    strikes, prices = [], []
    for i in range(n_strikes):
        x = -width * sd + 2.0 * width * sd * i / (n_strikes - 1)
        K = F * math.exp(x)
        v = vol_fn(K)
        # OTM option: put below the forward, call above.
        if K < F:
            prices.append(put_price(S0, K, t, r, v, b=r - q))
        else:
            prices.append(call_price(S0, K, t, r, v, b=r - q))
        strikes.append(K)
    return vix_from_chain(strikes, prices, F, t, r)


def svix_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=201, width=6.0):
    """Martin (2013) "simple variance" index (SVIX) from a smile ``vol_fn(K)``.

    Unlike the VIX log-contract (``1/K^2`` weights), the simple variance swap
    weights the OTM strip by ``1/F^2`` -- a constant -- so it corresponds to the
    payoff ``(S_T - F)^2 / F^2`` and needs no log approximation, making it robust
    to large moves/jumps and giving a genuine lower bound on the equity premium
    (Martin). The fair simple variance is

        SVIX^2 = (2 e^{r t} / (t F^2)) * ( OTM option strip ),

    reported as ``100 * SVIX``. A flat smile returns approximately
    ``100 * sigma`` (equal to VIX only to leading order; the two differ at higher
    order in vol).
    """
    from .bsm import call_price, put_price

    if t <= 0:
        raise ValueError("t must be positive")
    F = S0 * math.exp((r - q) * t)
    atm_vol = vol_fn(F)
    sd = atm_vol * math.sqrt(t)
    growth = math.exp(r * t)

    strikes, opts = [], []
    for i in range(n_strikes):
        x = -width * sd + 2.0 * width * sd * i / (n_strikes - 1)
        K = F * math.exp(x)
        v = vol_fn(K)
        if K < F:
            opts.append(put_price(S0, K, t, r, v, b=r - q))
        else:
            opts.append(call_price(S0, K, t, r, v, b=r - q))
        strikes.append(K)

    strip = 0.0
    for i in range(len(strikes) - 1):
        k0, k1 = strikes[i], strikes[i + 1]
        strip += 0.5 * (opts[i] + opts[i + 1]) * (k1 - k0)   # unweighted (1/F^2 outside)

    var = (2.0 * growth / (t * F * F)) * strip
    return var, 100.0 * math.sqrt(max(var, 0.0))


def equity_premium_lower_bound(S0, t, r, vol_fn, q=0.0, n_strikes=201,
                               width=6.0):
    """Martin's (2013) lower bound on the expected equity excess return.

    Martin shows that, under the (empirically mild) negative-correlation
    condition, the expected simple excess return of the market over ``[0, t]`` is
    bounded below by the risk-neutral *simple variance*:

        (1/t) E_0[ (R_market - R_f) ] >= Rf * SVIX^2,

    where ``SVIX^2`` is the annualized simple-variance index
    (:func:`svix_from_smile`) and ``Rf = e^{r t}`` the gross risk-free return.
    This returns the annualized lower bound ``Rf * SVIX^2`` -- a model-free floor
    on the equity premium computable purely from option prices.
    """
    var, _svix = svix_from_smile(S0, t, r, vol_fn, q=q, n_strikes=n_strikes,
                                 width=width)
    rf = math.exp(r * t)
    return rf * var
