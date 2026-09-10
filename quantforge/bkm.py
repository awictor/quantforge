"""Bakshi-Kapadia-Madan (2003) risk-neutral moments from an option smile.

BKM extract the risk-neutral variance, skewness and kurtosis of the log-return
model-free from a strip of out-of-the-money options. Three "contracts" -- the
quadratic (log), cubic and quartic payoffs -- are each spanned by an option
integral with a strike-dependent weight:

    V = integral w_V(K) O(K) dK,   W = integral w_W(K) O(K) dK,
    X = integral w_X(K) O(K) dK,

with (for K > S0, using calls; K < S0, using puts)

    w_V = 2(1 - ln(K/S0)) / K^2,
    w_W = (6 ln(K/S0) - 3 ln(K/S0)^2) / K^2,
    w_X = (12 ln(K/S0)^2 - 4 ln(K/S0)^3) / K^2.

The de-meaned risk-neutral moments then follow (Bakshi-Kapadia-Madan eqs. 5-8):

    mu = e^{rt} - 1 - e^{rt}(V/2 + W/6 + X/24),
    Var = e^{rt} V - mu^2,
    Skew = (e^{rt} W - 3 mu e^{rt} V + 2 mu^3) / Var^{1.5},
    Kurt = (e^{rt} X - 4 mu e^{rt} W + 6 mu^2 e^{rt} V - 3 mu^4) / Var^2.

A symmetric smile gives zero skewness; a downward (equity) skew gives negative
skewness and fat-tailed excess kurtosis. Pure standard library.
"""

import math


def bkm_moments_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401, width=8.0):
    """Risk-neutral (variance, skewness, excess kurtosis) via Bakshi-Kapadia-Madan.

    ``vol_fn(K)`` is the implied-vol smile; OTM options are priced with
    Black-Scholes and the three moment contracts integrated by the trapezoidal
    rule. Returns ``(variance, skewness, excess_kurtosis)`` of the ``t``-horizon
    log-return under the risk-neutral measure.
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
        # OTM relative to spot S0 (BKM references log(K/S0)).
        if K < S0:
            opts.append(put_price(S0, K, t, r, v, b=r - q))
        else:
            opts.append(call_price(S0, K, t, r, v, b=r - q))
        strikes.append(K)

    def w_V(K):
        lk = math.log(K / S0)
        return 2.0 * (1.0 - lk) / (K * K)

    def w_W(K):
        lk = math.log(K / S0)
        return (6.0 * lk - 3.0 * lk * lk) / (K * K)

    def w_X(K):
        lk = math.log(K / S0)
        return (12.0 * lk * lk - 4.0 * lk ** 3) / (K * K)

    def integ(weight):
        total = 0.0
        for i in range(len(strikes) - 1):
            k0, k1 = strikes[i], strikes[i + 1]
            g0 = weight(k0) * opts[i]
            g1 = weight(k1) * opts[i + 1]
            total += 0.5 * (g0 + g1) * (k1 - k0)
        return total

    V = integ(w_V)
    W = integ(w_W)
    X = integ(w_X)

    mu = growth - 1.0 - growth * (V / 2.0 + W / 6.0 + X / 24.0)
    var = growth * V - mu * mu
    if var <= 0.0:
        return var, 0.0, 0.0
    skew = (growth * W - 3.0 * mu * growth * V + 2.0 * mu ** 3) / var ** 1.5
    kurt = (growth * X - 4.0 * mu * growth * W
            + 6.0 * mu * mu * growth * V - 3.0 * mu ** 4) / (var * var)
    return var, skew, kurt - 3.0


def skew_swap_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401, width=8.0):
    """Fair skew-swap value: the risk-neutral skewness from the BKM moments."""
    _var, skew, _ek = bkm_moments_from_smile(S0, t, r, vol_fn, q, n_strikes,
                                             width)
    return skew
