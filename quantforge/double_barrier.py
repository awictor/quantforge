"""Double-barrier knock-out options (Kunitomo-Ikeda 1992), closed form.

A double knock-out call pays the vanilla call payoff at expiry only if the spot
stays strictly inside a corridor ``(L, U)`` for the whole life; touching either
the lower barrier ``L`` or the upper barrier ``U`` extinguishes it. Kunitomo and
Ikeda (1992) price it as a rapidly converging series of image terms that enforce
the two absorbing boundaries:

    C = sum_{n=-inf}^{inf} [ mu1_n * (BSM-type call term) - mu2_n * (...) ]

Here we use the standard formulation with flat barriers (no curvature) and carry
``b``. The series converges geometrically; a handful of terms is machine-accurate.

Pure standard library.
"""

import math

from .mathfns import norm_cdf


def double_knockout_call(S, K, L, U, t, r, sigma, b=None, q=0.0, terms=8):
    """Price a double-barrier knock-out call in closed form (Kunitomo-Ikeda).

    Parameters
    ----------
    S, K : float
        Spot and strike.
    L, U : float
        Lower and upper knock-out barriers with ``L < S < U``. Both must be
        positive and ``L < U``.
    t, r, sigma : float
        Time to expiry (years), risk-free rate, volatility.
    b : float, optional
        Cost of carry. Defaults to ``r - q``.
    q : float
        Continuous dividend yield, used only when ``b`` is not given.
    terms : int
        Number of image terms on each side of the series (total ``2*terms+1``).

    Returns
    -------
    float
        Value of the down-and-out-and-up-and-out call. Non-negative and never
        exceeds the vanilla call; it approaches the vanilla as the barriers move
        far away.
    """
    if S <= 0 or K <= 0 or t <= 0 or sigma <= 0:
        raise ValueError("S, K, t, sigma must be positive")
    if not (0.0 < L < S < U):
        raise ValueError("require 0 < L < S < U")
    if b is None:
        b = r

    vt = sigma * math.sqrt(t)
    var = sigma * sigma
    F = U / L
    # Kunitomo-Ikeda with flat barriers: delta1 = delta2 = 0. Only the curvature
    # parameters vanish; the image series over the log-corridor width remains.
    mu = b - 0.5 * var

    def d(x):
        return (math.log(x) + (b + 0.5 * var) * t) / vt

    total = 0.0
    for n in range(-terms, terms + 1):
        Fn = F ** n
        # Exponents for the two image reflections.
        d1 = (math.log(S * Fn * Fn / K) + (b + 0.5 * var) * t) / vt
        d2 = (math.log(S * Fn * Fn / U) + (b + 0.5 * var) * t) / vt
        d3 = (math.log(L * L / (S * K) * Fn * Fn) + (b + 0.5 * var) * t) / vt
        d4 = (math.log(L * L / (S * U) * Fn * Fn) + (b + 0.5 * var) * t) / vt

        pow1 = (Fn) ** (2.0 * (mu / var + 1.0))
        pow2 = (L / S * Fn) ** (2.0 * (mu / var + 1.0))
        # Value-side coefficients.
        c1 = (Fn) ** (2.0 * mu / var)
        c2 = (L / S * Fn) ** (2.0 * mu / var)

        term_call = (
            S * math.exp((b - r) * t) * (pow1 * (norm_cdf(d1) - norm_cdf(d2))
                                         - pow2 * (norm_cdf(d3) - norm_cdf(d4)))
            - K * math.exp(-r * t) * (c1 * (norm_cdf(d1 - vt) - norm_cdf(d2 - vt))
                                      - c2 * (norm_cdf(d3 - vt) - norm_cdf(d4 - vt)))
        )
        total += term_call
    return max(total, 0.0)
