"""CMS convexity adjustment: standard-model formula and static replication.

A constant-maturity-swap (CMS) payment pays a swap rate at a date other than the
swap's natural (annuity-weighted) settlement, so its expectation under the
payment measure is *not* the forward swap rate -- it carries a convexity
adjustment. Two routes are provided:

  * ``cms_adjustment_standard``: the classic linear-TSR / Hagan "standard model"
    closed form, ``CA = S0^2 * sigma^2 * T * G``, where ``sigma`` is the (Black)
    swap-rate vol, ``T`` the fixing time and ``G`` a level-function factor
    depending on the swap tenor and payment lag.
  * ``cms_rate_replication``: static replication of the CMS payoff by a strip of
    payer and receiver swaptions across strikes, using a caller-supplied
    volatility smile (e.g. from :class:`quantforge.VolCube`). This captures the
    smile, and with a flat smile it reproduces the standard-model adjustment.

Pure standard library.
"""

import math

from .bsm import call_price, put_price, OptionType


def _level_G(forward, tenor, freq, pay_lag):
    """Hagan linear-TSR level factor G for the standard CMS adjustment.

    Uses the standard bond-math weight of a par swap: with ``m = freq * tenor``
    payments the annuity derivative gives
    ``G = tenor / (1 - (1+S/ freq)^{-m})  -  ... ``; we use the widely-quoted
    closed form for an in-arrears/lagged CMS (Hagan 2003, linear TSR).
    """
    m = freq * tenor
    x = 1.0 + forward / freq
    # Annuity A(S) = (1 - x^{-m}) / S ; the convexity factor is
    # G = -A''(S)/A'(S) evaluated at the forward, plus the payment-lag delay.
    # Closed form of the level-function log-derivative (Hagan standard model):
    dA = (m / freq) * x ** (-m - 1) / forward \
        - (1.0 - x ** (-m)) / (forward * forward)
    A = (1.0 - x ** (-m)) / forward
    g = -dA / A
    # Payment-lag discounting delay term ~ pay_lag * 1/(1+S*pay_lag).
    g += pay_lag / (1.0 + forward * pay_lag)
    return g


def cms_adjustment_standard(forward, sigma, expiry, tenor, freq=1.0, pay_lag=0.0):
    """Standard-model (linear-TSR) CMS convexity adjustment.

    Args:
        forward: forward swap rate S0.
        sigma: lognormal (Black) swap-rate volatility.
        expiry: fixing time in years.
        tenor: swap tenor in years.
        freq: payment frequency of the underlying swap (per year).
        pay_lag: payment delay in years (0 = natural payment).

    Under the linear terminal-swap-rate model the adjustment is
    ``CA = G * Var_A(S_T) = G * S0^2 (e^{sigma^2 T} - 1)`` -- the exact
    lognormal variance, not just its ``sigma^2 T`` leading term. Returns the
    additive adjustment so ``E_pay[S_T] = forward + CA``.
    """
    G = _level_G(forward, tenor, freq, pay_lag)
    var = forward * forward * (math.exp(sigma * sigma * expiry) - 1.0)
    return G * var


def cms_rate_convexity_replication(forward, expiry, tenor, vol_fn, freq=1.0,
                                   pay_lag=0.0, width=8.0, n=800):
    """CMS convexity adjustment by static replication over a swaption strip.

    Under the linear-TSR model the adjustment is ``G * E_A[(S_T - S0)^2]``, and
    the second moment of the swap rate is replicated model-free by a strip of
    swaptions (Carr-Madan variance replication):

        E_A[(S_T - S0)^2] = 2 * integral_0^inf swaption(K) dK,

    with payer swaptions for ``K >= S0`` and receiver swaptions for ``K < S0``,
    each priced at the smile vol ``vol_fn(K)``. This captures the whole smile,
    and with a flat ``vol_fn`` it reproduces :func:`cms_adjustment_standard`.
    """
    G = _level_G(forward, tenor, freq, pay_lag)
    atm_vol = vol_fn(forward)
    sd = atm_vol * math.sqrt(expiry) * forward
    lo = max(forward - width * sd, 1e-8)
    hi = forward + width * sd
    dK = (hi - lo) / n

    def swaption(K):
        v = vol_fn(K)
        if K >= forward:
            return call_price(forward, K, expiry, 0.0, v, b=0.0)
        return put_price(forward, K, expiry, 0.0, v, b=0.0)

    total = 0.0
    prev = None
    for i in range(n + 1):
        Kc = lo + i * dK
        val = swaption(Kc)
        if prev is not None:
            total += 0.5 * (prev + val) * dK
        prev = val
    second_moment = 2.0 * total
    return G * second_moment


def cms_rate(forward, sigma, expiry, tenor, freq=1.0, pay_lag=0.0):
    """Convexity-adjusted expected CMS rate under the standard model."""
    return forward + cms_adjustment_standard(forward, sigma, expiry, tenor,
                                             freq, pay_lag)


def cms_adjustment_greeks(forward, sigma, expiry, tenor, freq=1.0, pay_lag=0.0):
    """Sensitivities of the standard-model CMS convexity adjustment.

    Central finite differences of :func:`cms_adjustment_standard` for
    ``d_forward`` (d(CA)/d forward) and ``d_sigma`` (d(CA)/d sigma). The
    adjustment is monotone increasing in the vol (more convexity), so
    ``d_sigma > 0``. Returns a dict with ``adjustment``, ``d_forward``,
    ``d_sigma``.
    """
    def ca(f=forward, sig=sigma):
        return cms_adjustment_standard(f, sig, expiry, tenor, freq, pay_lag)

    base = ca()
    hf = 1e-6 * max(abs(forward), 1e-4)
    d_forward = (ca(f=forward + hf) - ca(f=forward - hf)) / (2.0 * hf)
    hv = 1e-6
    d_sigma = (ca(sig=sigma + hv) - ca(sig=max(sigma - hv, 0.0))) / (
        (2.0 * hv) if sigma - hv >= 0 else hv)
    return {"adjustment": base, "d_forward": d_forward, "d_sigma": d_sigma}
