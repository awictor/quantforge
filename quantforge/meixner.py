"""Meixner Levy process option pricing (Schoutens 2002).

The Meixner process is a pure-jump Levy process built from the Meixner
distribution, a flexible three-parameter family that fits equity return skew and
excess kurtosis at least as well as NIG while keeping a fully analytic
characteristic function:

    phi(u) = ( cos(b/2) / cosh((a u - i b) / 2) )^{2 d}

so the characteristic *exponent* (with ``E[e^{i u L_t}] = e^{t psi(u)}``) is

    psi(u) = 2 d * ( ln cos(b/2) - ln cosh((a u - i b) / 2) ).

``a > 0`` scales the jump sizes, ``b in (-pi, pi)`` sets the asymmetry
(``b < 0`` gives the equity down-skew), and ``d > 0`` scales the activity. The
mean-correcting martingale drift ``omega = -psi(-i)`` is supplied by the shared
Carr-Madan engine, so Meixner is priced (and cross-checked by the COS method)
through exactly the same path as CGMY and NIG. Pure standard library.
"""

import cmath
import math

from .bsm import OptionType
from .carrmadan import levy_price


def _meixner_psi(u, a, b, d):
    """Meixner characteristic exponent psi(u) (u complex)."""
    return 2.0 * d * (cmath.log(math.cos(b / 2.0))
                      - cmath.log(cmath.cosh((a * u - 1j * b) / 2.0)))


def meixner_price(S, K, t, r, a, b, d, option_type=OptionType.CALL, q=0.0,
                  cm_alpha=1.5, upper=200.0) -> float:
    """Price a European option under the Meixner model via Carr-Madan inversion.

    Args:
        a: jump-size scale (> 0).
        b: asymmetry in ``(-pi, pi)``; ``b < 0`` gives a downward skew.
        d: activity / tail parameter (> 0).
        q: continuous dividend yield.
        cm_alpha: Carr-Madan damping. The transform needs the moment-generating
            function at ``s = cm_alpha + 1`` to be finite, i.e.
            ``a (cm_alpha + 1) + b < pi`` (the ``cosh`` argument must stay off its
            pole).

    Puts use put-call parity.
    """
    if a <= 0 or d <= 0:
        raise ValueError("a and d must be positive")
    if not (-math.pi < b < math.pi):
        raise ValueError("b must be in (-pi, pi)")
    # The MGF exists for real s with |a s + b| < pi; enforce it at s = cm_alpha+1
    # (and its reflection) so the damped transform converges.
    s = cm_alpha + 1.0
    if abs(a * s + b) >= math.pi:
        raise ValueError("need |a (cm_alpha + 1) + b| < pi for a finite "
                         "transform; lower cm_alpha or a")

    return levy_price(S, K, t, r, q,
                      lambda u: _meixner_psi(u, a, b, d),
                      option_type, alpha=cm_alpha, upper=upper)


def meixner_greeks(S, K, t, r, a, b, d, option_type=OptionType.CALL, q=0.0,
                   cm_alpha=1.5):
    """Greeks of a Meixner option by central finite differences.

    Central differences of :func:`meixner_price` for the spot Greeks ``delta``
    (dV/dS), ``gamma`` (d2V/dS2), and ``theta`` (calendar decay), plus the
    asymmetry (skew) sensitivity ``d_b`` (dV/db). The ``d_b`` bump is clipped to
    keep ``b`` in ``(-pi, pi)`` on both sides. Returns a dict with ``price``,
    ``delta``, ``gamma``, ``theta``, ``d_b``.
    """
    if a <= 0 or d <= 0:
        raise ValueError("a and d must be positive")
    if not (-math.pi < b < math.pi):
        raise ValueError("b must be in (-pi, pi)")
    if t <= 0:
        raise ValueError("t must be positive")

    def px(S_=S, t_=t, b_=b):
        return meixner_price(S_, K, t_, r, a, b_, d, option_type, q=q,
                             cm_alpha=cm_alpha)

    base = px()
    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    hb = min(1e-3, 0.25 * (math.pi - abs(b)))
    d_b = (px(b_=b + hb) - px(b_=b - hb)) / (2.0 * hb)
    return {"price": base, "delta": delta, "gamma": gamma, "theta": theta,
            "d_b": d_b}


def meixner_smile(S, strikes, t, r, a, b, d, q=0.0, cm_alpha=1.5):
    """Black-Scholes implied-vol smile the Meixner model produces.

    Prices a call at each strike and inverts to a Black-Scholes implied vol,
    returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
    ``F = S e^{(r-q) t}``. ``b < 0`` tilts the smile into a downward skew.
    """
    from .implied import implied_volatility

    F = S * math.exp((r - q) * t)
    out = []
    for K in sorted(strikes):
        c = meixner_price(S, K, t, r, a, b, d, OptionType.CALL, q=q,
                          cm_alpha=cm_alpha)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
