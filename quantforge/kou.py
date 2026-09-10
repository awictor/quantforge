"""Kou (2002) double-exponential jump-diffusion option pricing.

Kou replaces Merton's Gaussian jumps with an asymmetric double-exponential jump
size, which fattens both tails while allowing the up- and down-jump tails to
differ -- capturing the empirical negative skew and leptokurtosis of equity
returns with a still-analytic characteristic function. The log-return jump ``Y``
has density

    f(y) = p * eta1 * e^{-eta1 y}   (y >= 0, up jumps)
         + (1-p) * eta2 * e^{eta2 y} (y <  0, down jumps)

with ``eta1 > 1`` (so ``E[e^Y]`` is finite), ``eta2 > 0`` and up-probability
``p in [0, 1]``. Under the risk-neutral measure the log-spot is a Levy process

    x_t = ln S + (r - q - sigma^2/2 - lambda*zeta) t + sigma W_t + sum_{i<=N_t} Y_i,

where ``N_t`` is Poisson(``lambda``) and ``zeta = E[e^Y] - 1`` compensates the
jump drift so the discounted spot is a martingale. Pricing uses the same
two-probability Gauss-Legendre Fourier integral as the Heston/Bates pricers,
applied to the Kou log-price characteristic function. Pure standard library.
"""

import cmath
import math

from .bsm import OptionType, _coerce_type
from .heston import _GL_NODES, _GL_WEIGHTS


def _kou_kappa(u, p, eta1, eta2):
    """Jump moment-generating factor ``E[e^{u Y}]`` for a Kou jump (u complex).

    ``= p * eta1/(eta1 - u) + (1-p) * eta2/(eta2 + u)``; valid for
    ``Re(u) < eta1`` and ``Re(u) > -eta2``.
    """
    return p * eta1 / (eta1 - u) + (1.0 - p) * eta2 / (eta2 + u)


def _kou_psi(u, r, q, sigma, lam, p, eta1, eta2, zeta):
    """Levy characteristic exponent ``psi(u)`` with ``cf(u) = e^{t psi(u)}``."""
    return (u * (r - q - 0.5 * sigma * sigma - lam * zeta)
            + 0.5 * sigma * sigma * u * u
            + lam * (_kou_kappa(u, p, eta1, eta2) - 1.0))


def _kou_cf(phi, x0, t, r, q, sigma, lam, p, eta1, eta2, zeta):
    """Characteristic function of the log-spot ``E[e^{i phi x_t}]``."""
    return cmath.exp(1j * phi * x0 + t * _kou_psi(1j * phi, r, q, sigma,
                                                  lam, p, eta1, eta2, zeta))


def _kou_probability(S, K, t, r, q, sigma, lam, p, eta1, eta2, zeta, j,
                     upper=200.0):
    """Exercise probability ``Pi_j`` via Gauss-Legendre Fourier inversion."""
    x0 = math.log(S)
    lnK = math.log(K)
    fwd_mgf = S * math.exp((r - q) * t)   # E[S_T] = cf(-i)
    half = 0.5 * upper
    total = 0.0
    for node, w in zip(_GL_NODES, _GL_WEIGHTS):
        phi = half * (node + 1.0)
        if phi <= 0:
            phi = 1e-8
        if j == 2:
            f = _kou_cf(phi, x0, t, r, q, sigma, lam, p, eta1, eta2, zeta)
        else:
            # Share-measure density: cf shifted by -i, normalised by the forward.
            f = (_kou_cf(phi - 1j, x0, t, r, q, sigma, lam, p, eta1, eta2, zeta)
                 / fwd_mgf)
        integrand = (cmath.exp(-1j * phi * lnK) * f / (1j * phi)).real
        total += w * integrand
    return 0.5 + half * total / math.pi


def kou_price(S, K, t, r, sigma, lam, p, eta1, eta2,
              option_type=OptionType.CALL, q=0.0, upper=200.0) -> float:
    """Price a European option under the Kou double-exponential jump-diffusion.

    Args:
        sigma: diffusion volatility.
        lam: jump intensity (expected jumps per year, >= 0).
        p: probability a jump is upward (in [0, 1]).
        eta1: up-jump tail rate (must be > 1 so E[e^Y] is finite).
        eta2: down-jump tail rate (> 0).
        q: continuous dividend yield.

    ``lam = 0`` recovers Black-Scholes. Puts use put-call parity.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if sigma < 0 or lam < 0:
        raise ValueError("sigma and lam must be non-negative")
    if not (0.0 <= p <= 1.0):
        raise ValueError("p must be in [0, 1]")
    if eta1 <= 1.0 or eta2 <= 0.0:
        raise ValueError("need eta1 > 1 and eta2 > 0")

    if t == 0:
        return max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)

    zeta = p * eta1 / (eta1 - 1.0) + (1.0 - p) * eta2 / (eta2 + 1.0) - 1.0

    P1 = _kou_probability(S, K, t, r, q, sigma, lam, p, eta1, eta2, zeta, 1, upper)
    P2 = _kou_probability(S, K, t, r, q, sigma, lam, p, eta1, eta2, zeta, 2, upper)

    call = S * math.exp(-q * t) * P1 - K * math.exp(-r * t) * P2
    if ot is OptionType.CALL:
        return call
    return call - S * math.exp(-q * t) + K * math.exp(-r * t)


def kou_smile(S, strikes, t, r, sigma, lam, p, eta1, eta2, q=0.0):
    """Black-Scholes implied-vol smile the Kou model produces.

    Prices a call at each strike and inverts to a Black-Scholes implied vol,
    returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
    ``F = S e^{(r-q) t}``. An asymmetric jump distribution (``eta1 != eta2`` or
    ``p != 1/2``) tilts the smile into a skew.
    """
    from .implied import implied_volatility

    F = S * math.exp((r - q) * t)
    out = []
    for K in sorted(strikes):
        c = kou_price(S, K, t, r, sigma, lam, p, eta1, eta2,
                      OptionType.CALL, q=q)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
