"""CGMY (Carr-Geman-Madan-Yor 2002) tempered-stable Levy option pricing.

CGMY is a pure-jump Levy process whose Levy density is a stable density with
independent exponential tempering on each side:

    nu(x) = C * e^{-G|x|} / |x|^{1+Y}   (x < 0, down jumps)
          = C * e^{-M x}   / x^{1+Y}    (x > 0, up jumps)

``C > 0`` scales overall jump activity, ``G, M > 0`` temper the down/up tails
(``G != M`` gives skew), and the fine-structure exponent ``Y < 2`` sets the path
behaviour: ``Y < 0`` finite activity, ``0 <= Y < 1`` infinite activity / finite
variation, ``1 <= Y < 2`` infinite variation. ``Y -> 0`` recovers Variance
Gamma. The characteristic exponent is closed form,

    psi(u) = C * Gamma(-Y) * [ (M - i u)^Y - M^Y + (G + i u)^Y - G^Y ],

so options are priced by Carr-Madan Fourier inversion of the damped call
transform, integrated here with the shared Gauss-Legendre nodes. Pure standard
library (uses ``cmath`` for the complex powers; ``Gamma(-Y)`` is real).
"""

import math

from .bsm import OptionType
from .carrmadan import levy_price


def _cgmy_psi(u, C, G, M, Y):
    """CGMY characteristic exponent psi(u): cf(u) = exp(t * psi(u)) (u complex)."""
    gamma_negY = math.gamma(-Y)   # real for non-integer Y in (0, 2)
    return C * gamma_negY * ((M - 1j * u) ** Y - M ** Y
                             + (G + 1j * u) ** Y - G ** Y)


def _cgmy_char_logspot(u, x0, t, r, q, C, G, M, Y):
    """Characteristic function of ln S_T under the risk-neutral CGMY measure.

    Retained for reference and independent (e.g. Gil-Pelaez) cross-checks; the
    Carr-Madan pricer itself now runs through :mod:`quantforge.carrmadan`.
    """
    import cmath
    omega = -_cgmy_psi(-1j, C, G, M, Y)
    drift = x0 + (r - q + omega) * t
    return cmath.exp(1j * u * drift + t * _cgmy_psi(u, C, G, M, Y))


def cgmy_price(S, K, t, r, C, G, M, Y, option_type=OptionType.CALL, q=0.0,
               alpha=1.5, upper=200.0) -> float:
    """Price a European option under the CGMY model by Carr-Madan inversion.

    Args:
        C, G, M, Y: CGMY parameters (``C > 0``; ``G, M > 0``; ``Y < 2`` and not a
            non-negative integer -- ``Y = 0`` is Variance Gamma, handled by the
            limit only approximately here, so pass a small ``Y`` instead).
        alpha: Carr-Madan damping factor (> 0); the call transform needs
            ``E[S_T^{alpha+1}] < infinity``, i.e. ``alpha + 1 < M``.
        upper: Fourier-integral truncation.

    ``C = 0`` gives a degenerate (deterministic-forward) payoff. Puts use
    put-call parity.
    """
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if C < 0 or G <= 0 or M <= 0:
        raise ValueError("need C >= 0, G > 0, M > 0")
    if Y >= 2.0:
        raise ValueError("Y must be < 2")
    if alpha <= 0 or alpha + 1.0 >= M:
        raise ValueError("need 0 < alpha and alpha + 1 < M for a finite transform")

    return levy_price(S, K, t, r, q,
                      lambda u: _cgmy_psi(u, C, G, M, Y),
                      option_type, alpha=alpha, upper=upper)


def cgmy_smile(S, strikes, t, r, C, G, M, Y, q=0.0, alpha=1.5):
    """Black-Scholes implied-vol smile the CGMY model produces.

    Prices a call at each strike and inverts to a Black-Scholes implied vol,
    returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
    ``F = S e^{(r-q) t}``. Tail asymmetry (``G != M``) tilts the smile into a
    skew; smaller ``Y`` fattens the wings.
    """
    from .implied import implied_volatility

    F = S * math.exp((r - q) * t)
    out = []
    for K in sorted(strikes):
        c = cgmy_price(S, K, t, r, C, G, M, Y, OptionType.CALL, q=q, alpha=alpha)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
