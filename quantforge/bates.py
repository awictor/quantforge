"""Bates (1996) model: Heston stochastic volatility plus Merton lognormal jumps.

The Bates model adds Merton's compound-Poisson lognormal jumps to the Heston
square-root variance process, so the log-spot has both a diffusive stochastic-
vol component and jump risk:

    dS/S = (r - q - lambda*k) dt + sqrt(V) dW_S + (e^J - 1) dN
    dV   = kappa (theta - V) dt + xi sqrt(V) dW_V,   d<W_S, W_V> = rho dt

with jumps arriving at Poisson rate ``lambda``, jump log-sizes
``J ~ N(mu_j, sigma_j^2)``, mean proportional jump ``k = e^{mu_j + sigma_j^2/2} - 1``,
and a drift compensator ``-lambda*k`` keeping the discounted spot a martingale.

Because the jumps are independent of the diffusion, the Bates characteristic
function is the Heston one multiplied by the Merton jump characteristic factor;
the option is priced by the same two-probability Fourier integral as Heston,
reusing its Gauss-Legendre machinery. Pure standard library.
"""

import cmath
import math

from .bsm import OptionType, _coerce_type
from .heston import _GL_NODES, _GL_WEIGHTS


def _jump_cf_exponent(phi, t, lam, mu_j, sigma_j, j):
    """Log of the Merton-jump characteristic factor for probability ``j``.

    In the Heston two-probability decomposition, ``P2`` integrates the log-spot
    characteristic function under the risk-neutral measure (argument ``phi``)
    while ``P1`` integrates it under the share measure, whose argument is shifted
    by ``-i`` (i.e. ``phi - i``). Writing ``a = i*phi + s`` with ``s = 1`` for
    ``j = 1`` (the ``-i`` shift: ``i(phi - i) = i*phi + 1``) and ``s = 0`` for
    ``j = 2``, the compound-Poisson lognormal jump contributes to
    ``ln E[e^{a X_T}]``

        lam*t*( E[e^{a J}] - 1 ) - a * lam * t * k ,

    where the last term is the martingale drift compensator and
    ``k = e^{mu_j + sigma_j^2/2} - 1`` is the mean proportional jump.
    """
    s = 1.0 if j == 1 else 0.0
    k = math.exp(mu_j + 0.5 * sigma_j * sigma_j) - 1.0
    a = 1j * phi + s
    # E[e^{a J}] for J ~ N(mu_j, sigma_j^2).
    ej = cmath.exp(a * mu_j + 0.5 * a * a * sigma_j * sigma_j)
    return lam * t * (ej - 1.0) - a * lam * t * k


def _bates_char(phi, S, K, t, r, q, v0, kappa, theta, xi, rho,
                lam, mu_j, sigma_j, j):
    """Bates characteristic integrand: Heston 'little trap' times the jump factor."""
    x = math.log(S)
    a = kappa * theta
    if j == 1:
        u = 0.5
        b = kappa - rho * xi
    else:
        u = -0.5
        b = kappa

    rxi = rho * xi
    d = cmath.sqrt((rxi * phi * 1j - b) ** 2 - xi * xi * (2.0 * u * phi * 1j - phi * phi))
    g = (b - rxi * phi * 1j - d) / (b - rxi * phi * 1j + d)
    exp_dt = cmath.exp(-d * t)

    C = ((r - q) * phi * 1j * t
         + (a / (xi * xi)) * ((b - rxi * phi * 1j - d) * t
                              - 2.0 * cmath.log((1.0 - g * exp_dt) / (1.0 - g))))
    D = ((b - rxi * phi * 1j - d) / (xi * xi)) * ((1.0 - exp_dt) / (1.0 - g * exp_dt))

    jump = _jump_cf_exponent(phi, t, lam, mu_j, sigma_j, j)
    f = cmath.exp(C + D * v0 + jump + 1j * phi * x)
    return (cmath.exp(-1j * phi * math.log(K)) * f / (1j * phi)).real


def _bates_probability(S, K, t, r, q, v0, kappa, theta, xi, rho,
                       lam, mu_j, sigma_j, j, upper=200.0):
    half = 0.5 * upper
    total = 0.0
    for node, w in zip(_GL_NODES, _GL_WEIGHTS):
        phi = half * (node + 1.0)
        if phi <= 0:
            phi = 1e-8
        total += w * _bates_char(phi, S, K, t, r, q, v0, kappa, theta, xi, rho,
                                 lam, mu_j, sigma_j, j)
    return 0.5 + half * total / math.pi


def bates_price(S, K, t, r, v0, kappa, theta, xi, rho,
                lam, mu_j, sigma_j, option_type=OptionType.CALL,
                q=0.0, upper=200.0) -> float:
    """Price a European option under the Bates (Heston + Merton jumps) model.

    Args:
        v0, kappa, theta, xi, rho: Heston stochastic-variance parameters.
        lam: jump intensity (expected jumps per year, >= 0).
        mu_j, sigma_j: mean and std of the log jump size.
        q: continuous dividend yield.
        upper: Fourier-integral truncation (raise for long maturities/large xi).

    ``lam = 0`` recovers the Heston price exactly. Puts use put-call parity.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if v0 < 0 or theta < 0 or xi < 0:
        raise ValueError("variance parameters must be non-negative")
    if lam < 0 or sigma_j < 0:
        raise ValueError("lam and sigma_j must be non-negative")

    if t == 0:
        return max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)

    P1 = _bates_probability(S, K, t, r, q, v0, kappa, theta, xi, rho,
                            lam, mu_j, sigma_j, 1, upper)
    P2 = _bates_probability(S, K, t, r, q, v0, kappa, theta, xi, rho,
                            lam, mu_j, sigma_j, 2, upper)

    call = S * math.exp(-q * t) * P1 - K * math.exp(-r * t) * P2
    if ot is OptionType.CALL:
        return call
    return call - S * math.exp(-q * t) + K * math.exp(-r * t)


def bates_greeks(S, K, t, r, v0, kappa, theta, xi, rho, lam, mu_j, sigma_j,
                 option_type=OptionType.CALL, q=0.0):
    """Greeks of a Bates (Heston + jumps) option by central finite differences.

    Central differences of :func:`bates_price` for the spot Greeks ``delta``
    (dV/dS) and ``gamma`` (d2V/dS2), the initial-variance sensitivity ``vega_v0``
    (dV/dv0 -- the stochastic-vol analogue of vega), and the jump-intensity
    sensitivity ``d_lambda`` (dV/dlam). At ``lam = 0`` the Greeks reduce to the
    Heston Greeks. Returns a dict with ``price``, ``delta``, ``gamma``,
    ``vega_v0``, ``d_lambda``.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    if v0 < 0 or theta < 0 or xi < 0:
        raise ValueError("variance parameters must be non-negative")
    if lam < 0 or sigma_j < 0:
        raise ValueError("lam and sigma_j must be non-negative")

    def px(S_=S, v0_=v0, lam_=lam):
        return bates_price(S_, K, t, r, v0_, kappa, theta, xi, rho, lam_, mu_j,
                           sigma_j, ot, q=q)

    base = px()
    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega_v0 = (px(v0_=v0 + hv) - px(v0_=max(v0 - hv, 0.0))) / (
        (2.0 * hv) if v0 - hv >= 0 else hv)
    hl = 1e-3
    d_lambda = (px(lam_=lam + hl) - px(lam_=max(lam - hl, 0.0))) / (
        (2.0 * hl) if lam - hl >= 0 else hl)
    return {"price": base, "delta": delta, "gamma": gamma, "vega_v0": vega_v0,
            "d_lambda": d_lambda}


def bates_smile(S, strikes, t, r, v0, kappa, theta, xi, rho,
                lam, mu_j, sigma_j, q=0.0):
    """Black-Scholes implied-vol smile the Bates model produces.

    Prices a call at each strike and inverts to a Black-Scholes implied vol,
    returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
    ``F = S e^{(r-q) t}``. Jumps steepen the short-dated skew beyond what the
    Heston diffusion alone can produce.
    """
    from .implied import implied_volatility

    F = S * math.exp((r - q) * t)
    out = []
    for K in sorted(strikes):
        c = bates_price(S, K, t, r, v0, kappa, theta, xi, rho,
                        lam, mu_j, sigma_j, OptionType.CALL, q=q)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
