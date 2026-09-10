"""Variance-Gamma (Madan-Carr-Chang 1998) option pricing.

The Variance-Gamma model replaces the Brownian diffusion with Brownian motion
evaluated at a random, gamma-distributed business time. It is a pure-jump
process with three parameters:

    sigma : volatility of the Brownian component
    nu    : variance rate of the gamma time change (controls kurtosis)
    theta : drift of the Brownian component (controls skew; theta < 0 -> the
            equity-style left skew)

As ``nu -> 0`` the model collapses to Black-Scholes. The characteristic exponent
is closed form,

    psi(u) = -(1/nu) * log(1 - i theta nu u + 0.5 sigma^2 nu u^2),

so VG is priced through the shared Carr-Madan engine (:mod:`quantforge.carrmadan`)
and gets the COS-method cross-check for free, like the other Levy models.
"""

import cmath
import math

from .bsm import OptionType
from .carrmadan import levy_price


def _vg_psi(u, sigma, nu, theta):
    """Variance-Gamma characteristic exponent psi(u) (u complex)."""
    return -(1.0 / nu) * cmath.log(
        1.0 - 1j * theta * nu * u + 0.5 * sigma * sigma * nu * u * u)


def variance_gamma_price(S, K, t, r, sigma, nu, theta,
                         option_type=OptionType.CALL, q=0.0, cm_alpha=1.5,
                         upper=200.0):
    """Price a European option under the Variance-Gamma model.

    Args:
        sigma: Brownian volatility. nu: gamma-time variance rate (> 0).
        theta: Brownian drift (skew; negative for an equity left skew).
        q: continuous dividend yield.
        cm_alpha: Carr-Madan damping; the transform needs
            ``1 - theta nu (cm_alpha+1) - 0.5 sigma^2 nu (cm_alpha+1)^2 > 0``.

    Puts follow from put-call parity. As ``nu -> 0`` the price approaches the
    Black-Scholes value.
    """
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if nu <= 0 or sigma <= 0:
        raise ValueError("nu and sigma must be positive")
    # Martingale time-change condition: 1 - theta nu - 0.5 sigma^2 nu > 0.
    if 1.0 - theta * nu - 0.5 * sigma * sigma * nu <= 0:
        raise ValueError("parameters violate the VG martingale condition")
    # Carr-Madan damping needs the MGF finite at s = cm_alpha + 1.
    s = cm_alpha + 1.0
    if 1.0 - theta * nu * s - 0.5 * sigma * sigma * nu * s * s <= 0:
        raise ValueError("need a smaller cm_alpha: MGF diverges at cm_alpha + 1")

    return levy_price(S, K, t, r, q,
                      lambda u: _vg_psi(u, sigma, nu, theta),
                      option_type, alpha=cm_alpha, upper=upper)


def variance_gamma_greeks(S, K, t, r, sigma, nu, theta,
                          option_type=OptionType.CALL, q=0.0, cm_alpha=1.5):
    """Greeks of a Variance-Gamma option by central finite differences.

    Central differences of :func:`variance_gamma_price` for ``delta`` (dV/dS),
    ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma, the Brownian-vol sensitivity), and
    ``theta_greek`` (calendar decay, ``-dV/dt``). As ``nu -> 0`` the Greeks
    approach the Black-Scholes Greeks. ``theta`` is the VG skew *parameter*; the
    calendar Greek is returned as ``theta_greek`` to avoid the name clash.
    Returns a dict with ``price``, ``delta``, ``gamma``, ``vega``, ``theta_greek``.
    """
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    if nu <= 0 or sigma <= 0:
        raise ValueError("nu and sigma must be positive")

    def px(S_=S, t_=t, sigma_=sigma):
        return variance_gamma_price(S_, K, t_, r, sigma_, nu, theta,
                                    option_type, q=q, cm_alpha=cm_alpha)

    base = px()
    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta_greek = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta_greek": theta_greek}


def variance_gamma_smile(S, strikes, t, r, sigma, nu, theta, q=0.0,
                         cm_alpha=1.5):
    """Black-Scholes implied-vol smile the Variance-Gamma model produces.

    Prices a call at each strike and inverts to a Black-Scholes implied vol,
    returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
    ``F = S e^{(r-q) t}``. ``theta < 0`` tilts the smile into a downward skew;
    larger ``nu`` fattens the wings.
    """
    from .implied import implied_volatility

    F = S * math.exp((r - q) * t)
    out = []
    for K in sorted(strikes):
        c = variance_gamma_price(S, K, t, r, sigma, nu, theta,
                                 OptionType.CALL, q=q, cm_alpha=cm_alpha)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
