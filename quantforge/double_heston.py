"""Double-Heston (two-factor stochastic-volatility) option pricing.

Christoffersen, Heston & Jacobs (2009) show a single-factor Heston cannot match
both the level and the term structure of the smile: one mean-reversion speed
ties the short- and long-dated skew together. Double Heston adds a second,
independent variance process,

    dS/S = (r - q) dt + sqrt(V1) dW1 + sqrt(V2) dW2,
    dV_k = kappa_k (theta_k - V_k) dt + xi_k sqrt(V_k) dZ_k,  d<W_k, Z_k> = rho_k,

so a fast-reverting factor governs the short end and a slow one the long end.
Because the two factors are independent, the log-spot characteristic function is
the product of the two single-factor Heston pieces (sharing the drift once), and
the option is priced by the same Gil-Pelaez two-probability Fourier integral as
:mod:`quantforge.heston`, with the shared Gauss-Legendre nodes. Setting one
factor's vol-of-vol and initial variance to zero recovers single-factor Heston.
Pure standard library.
"""

import cmath
import math

from .bsm import OptionType, _coerce_type
from .heston import _GL_NODES, _GL_WEIGHTS


def _factor_exponent(u, t, v0, kappa, theta, xi, rho):
    """One Heston factor's contribution ``A + B v0`` to ln E[e^{i u X_t}]."""
    xi2 = xi * xi
    if xi2 < 1e-30:
        # Deterministic-variance limit: V_t = theta + (v0 - theta) e^{-kappa t}.
        # Its contribution is -0.5 (i u + u^2) * integral_0^t V_s ds.
        if abs(kappa) < 1e-12:
            iv = v0 * t
        else:
            iv = theta * t + (v0 - theta) * (1.0 - math.exp(-kappa * t)) / kappa
        return -0.5 * (1j * u + u * u) * iv
    d = cmath.sqrt((rho * xi * 1j * u - kappa) ** 2 + xi2 * (1j * u + u * u))
    num = kappa - rho * xi * 1j * u - d
    g = num / (kappa - rho * xi * 1j * u + d)
    exp_dt = cmath.exp(-d * t)
    B = (num / xi2) * ((1.0 - exp_dt) / (1.0 - g * exp_dt))
    A = (kappa * theta / xi2) * (num * t
                                 - 2.0 * cmath.log((1.0 - g * exp_dt) / (1.0 - g)))
    return A + B * v0


def _dh_cf(u, S, t, r, q, f1, f2):
    """Characteristic function of ln S_T under double Heston (u complex).

    ``f1``/``f2`` are ``(v0, kappa, theta, xi, rho)`` tuples for the two factors.
    """
    x = math.log(S)
    drift = 1j * u * (x + (r - q) * t)
    e1 = _factor_exponent(u, t, *f1)
    e2 = _factor_exponent(u, t, *f2)
    return cmath.exp(drift + e1 + e2)


def _dh_probability(S, K, t, r, q, f1, f2, j, upper=200.0):
    """Exercise probability ``Pi_j`` via Gil-Pelaez Fourier inversion."""
    lnK = math.log(K)
    fwd = _dh_cf(-1j, S, t, r, q, f1, f2)   # = S e^{(r-q)t}
    half = 0.5 * upper
    total = 0.0
    for node, w in zip(_GL_NODES, _GL_WEIGHTS):
        phi = half * (node + 1.0)
        if phi <= 0:
            phi = 1e-8
        if j == 1:
            cf = _dh_cf(phi - 1j, S, t, r, q, f1, f2) / fwd
        else:
            cf = _dh_cf(phi, S, t, r, q, f1, f2)
        total += w * (cmath.exp(-1j * phi * lnK) * cf / (1j * phi)).real
    return 0.5 + half * total / math.pi


def double_heston_price(S, K, t, r,
                        v01, kappa1, theta1, xi1, rho1,
                        v02, kappa2, theta2, xi2, rho2,
                        option_type=OptionType.CALL, q=0.0, upper=200.0) -> float:
    """Price a European option under the double-Heston model.

    Factor 1 is ``(v01, kappa1, theta1, xi1, rho1)`` and factor 2
    ``(v02, kappa2, theta2, xi2, rho2)`` -- typically a fast- and a slow-reverting
    variance. ``q`` is the dividend yield. Zeroing the second factor's ``xi2``
    and ``v02`` recovers single-factor Heston. Puts use put-call parity.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    for v0, th, xi in ((v01, theta1, xi1), (v02, theta2, xi2)):
        if v0 < 0 or th < 0 or xi < 0:
            raise ValueError("variance parameters must be non-negative")
    if t == 0:
        return max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)

    f1 = (v01, kappa1, theta1, xi1, rho1)
    f2 = (v02, kappa2, theta2, xi2, rho2)
    P1 = _dh_probability(S, K, t, r, q, f1, f2, 1, upper)
    P2 = _dh_probability(S, K, t, r, q, f1, f2, 2, upper)
    call = S * math.exp(-q * t) * P1 - K * math.exp(-r * t) * P2
    if ot is OptionType.CALL:
        return call
    return call - S * math.exp(-q * t) + K * math.exp(-r * t)


def double_heston_greeks(S, K, t, r,
                         v01, kappa1, theta1, xi1, rho1,
                         v02, kappa2, theta2, xi2, rho2,
                         option_type=OptionType.CALL, q=0.0):
    """Greeks of a double-Heston option by central finite differences.

    Central differences of :func:`double_heston_price` for the spot Greeks
    ``delta`` (dV/dS) and ``gamma`` (d2V/dS2), plus a per-factor
    initial-variance sensitivity ``vega_v01`` and ``vega_v02`` (dV/dv0 for each
    variance factor -- the stochastic-vol analogue of vega). Returns a dict with
    ``price``, ``delta``, ``gamma``, ``vega_v01``, ``vega_v02``.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    for v0, th, xi in ((v01, theta1, xi1), (v02, theta2, xi2)):
        if v0 < 0 or th < 0 or xi < 0:
            raise ValueError("variance parameters must be non-negative")

    def px(S_=S, dv1=0.0, dv2=0.0):
        return double_heston_price(S_, K, t, r,
                                   v01 + dv1, kappa1, theta1, xi1, rho1,
                                   v02 + dv2, kappa2, theta2, xi2, rho2,
                                   ot, q=q)

    base = px()
    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega_v01 = (px(dv1=hv) - px(dv1=-hv)) / (2.0 * hv) if v01 - hv >= 0 \
        else (px(dv1=hv) - base) / hv
    vega_v02 = (px(dv2=hv) - px(dv2=-hv)) / (2.0 * hv) if v02 - hv >= 0 \
        else (px(dv2=hv) - base) / hv
    return {"price": base, "delta": delta, "gamma": gamma,
            "vega_v01": vega_v01, "vega_v02": vega_v02}


def double_heston_smile(S, strikes, t, r,
                        v01, kappa1, theta1, xi1, rho1,
                        v02, kappa2, theta2, xi2, rho2, q=0.0):
    """Black-Scholes implied-vol smile the double-Heston model produces.

    Returns ``(log_moneyness, vol)`` pairs sorted by strike on the forward
    ``F = S e^{(r-q) t}``. The two mean-reversion speeds let the short- and
    long-dated skew move more independently than single-factor Heston allows.
    """
    from .implied import implied_volatility

    F = S * math.exp((r - q) * t)
    out = []
    for K in sorted(strikes):
        c = double_heston_price(S, K, t, r, v01, kappa1, theta1, xi1, rho1,
                                v02, kappa2, theta2, xi2, rho2,
                                OptionType.CALL, q=q)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
