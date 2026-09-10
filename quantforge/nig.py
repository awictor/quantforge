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
