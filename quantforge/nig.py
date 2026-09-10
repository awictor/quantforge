"""Normal Inverse Gaussian (NIG) Levy option pricing (Barndorff-Nielsen 1997).

NIG is a pure-jump Levy process obtained by subordinating Brownian motion to an
inverse-Gaussian clock. It is a flexible four-parameter family that fits equity
return skew and heavy tails with a simple closed-form characteristic exponent:

    psi(u) = delta * ( sqrt(alpha^2 - beta^2) - sqrt(alpha^2 - (beta + i u)^2) )

with tail heaviness ``alpha > 0``, asymmetry ``|beta| < alpha`` (``beta < 0``
gives the equity down-skew), and scale ``delta > 0``. Unlike CGMY, NIG has
infinite activity but a genuinely analytic exponent with no gamma function.

Priced by the shared Carr-Madan Fourier engine (:mod:`quantforge.carrmadan`),
which supplies the martingale drift correction ``omega = -psi(-i)`` and the
damped call inversion. Pure standard library.
"""

import cmath
import math

from .bsm import OptionType
from .carrmadan import levy_price


def _nig_psi(u, alpha, beta, delta):
    """NIG characteristic exponent psi(u) (u complex)."""
    a2 = alpha * alpha
    return delta * (cmath.sqrt(a2 - beta * beta)
                    - cmath.sqrt(a2 - (beta + 1j * u) ** 2))


def nig_price(S, K, t, r, alpha, beta, delta, option_type=OptionType.CALL,
              q=0.0, cm_alpha=1.5, upper=200.0) -> float:
    """Price a European option under the NIG model via Carr-Madan inversion.

    Args:
        alpha: tail-heaviness / steepness (> 0); larger = lighter tails.
        beta: asymmetry (``|beta| < alpha``); ``beta < 0`` gives a downward skew.
        delta: scale (> 0).
        q: continuous dividend yield.
        cm_alpha: Carr-Madan damping; needs ``alpha - (beta + cm_alpha + 1) > 0``
            for the martingale transform to stay finite.

    Puts use put-call parity.
    """
    if alpha <= 0 or delta <= 0:
        raise ValueError("alpha and delta must be positive")
    if abs(beta) >= alpha:
        raise ValueError("need |beta| < alpha")
    # The MGF at s = cm_alpha + 1 requires alpha^2 - (beta + s)^2 >= 0.
    s = cm_alpha + 1.0
    if alpha * alpha - (beta + s) ** 2 <= 0.0:
        raise ValueError("need alpha^2 > (beta + cm_alpha + 1)^2 for a finite "
                         "transform; lower cm_alpha or raise alpha")

    return levy_price(S, K, t, r, q,
                      lambda u: _nig_psi(u, alpha, beta, delta),
                      option_type, alpha=cm_alpha, upper=upper)


def nig_greeks(S, K, t, r, alpha, beta, delta, option_type=OptionType.CALL,
               q=0.0, cm_alpha=1.5):
    """Greeks of a NIG option by central finite differences.

    Central differences of :func:`nig_price` for the spot Greeks ``delta``
    (dV/dS), ``gamma`` (d2V/dS2), and ``theta`` (calendar decay), plus the
    process-parameter sensitivities ``d_alpha`` (dV/dalpha, tail steepness) and
    ``d_beta`` (dV/dbeta, skew). Returns a dict with ``price``, ``delta``,
    ``gamma``, ``theta``, ``d_alpha``, ``d_beta``.
    """
    if alpha <= 0 or delta <= 0:
        raise ValueError("alpha and delta must be positive")
    if abs(beta) >= alpha:
        raise ValueError("need |beta| < alpha")
    if t <= 0:
        raise ValueError("t must be positive")

    def px(S_=S, t_=t, a=alpha, bta=beta):
        return nig_price(S_, K, t_, r, a, bta, delta, option_type, q=q,
                         cm_alpha=cm_alpha)

    base = px()
    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    d_delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    ha = 1e-3 * alpha
    d_alpha = (px(a=alpha + ha) - px(a=alpha - ha)) / (2.0 * ha)
    hb = 1e-3 * max(abs(beta), 1.0)
    # Keep |beta| < alpha on both sides of the bump.
    hb = min(hb, 0.5 * (alpha - abs(beta)))
    d_beta = (px(bta=beta + hb) - px(bta=beta - hb)) / (2.0 * hb)
    return {"price": base, "delta": d_delta, "gamma": gamma, "theta": theta,
            "d_alpha": d_alpha, "d_beta": d_beta}


def nig_smile(S, strikes, t, r, alpha, beta, delta, q=0.0, cm_alpha=1.5):
    """Black-Scholes implied-vol smile the NIG model produces.

    Prices a call at each strike and inverts to a Black-Scholes implied vol,
    returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
    ``F = S e^{(r-q) t}``. ``beta < 0`` tilts the smile into a downward skew.
    """
    from .implied import implied_volatility

    F = S * math.exp((r - q) * t)
    out = []
    for K in sorted(strikes):
        c = nig_price(S, K, t, r, alpha, beta, delta, OptionType.CALL, q=q,
                      cm_alpha=cm_alpha)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
