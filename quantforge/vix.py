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
