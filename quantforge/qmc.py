"""Low-discrepancy (quasi-Monte Carlo) sequences and a QMC European pricer.

Pseudo-random Monte Carlo has O(1/sqrt(N)) error. Low-discrepancy sequences
fill the unit cube more evenly and give roughly O((log N)^d / N) error, so a
QMC estimate converges markedly faster for low-dimensional problems like a
single-step European payoff.

We use the Halton sequence (van der Corput radical inverse in a distinct prime
base per dimension) — exact, deterministic, and dependency-free. Uniforms are
mapped to normals with the accurate inverse-CDF from :mod:`quantforge.mathfns`,
so the pricer is a genuine QMC integration of the Black-Scholes payoff (not a
random simulation).
"""

import math
from typing import List

from .mathfns import norm_ppf
from .bsm import OptionType, _coerce_type, _validate, price as bsm_price

_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def _radical_inverse(n: int, base: int) -> float:
    """Van der Corput radical inverse of ``n`` in the given base, in [0, 1)."""
    inv = 0.0
    f = 1.0 / base
    while n > 0:
        inv += (n % base) * f
        n //= base
        f /= base
    return inv


def halton(index: int, dim: int) -> List[float]:
    """The ``index``-th Halton point in ``dim`` dimensions (0-based index).

    Skips index 0 (the origin) by convention via a 1-based offset internally,
    so callers can pass 0, 1, 2, ... and get well-spread points.
    """
    if dim < 1 or dim > len(_PRIMES):
        raise ValueError(f"dim must be in 1..{len(_PRIMES)}")
    n = index + 1
    return [_radical_inverse(n, _PRIMES[d]) for d in range(dim)]


def european_qmc(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                 n_points=8192) -> float:
    """Price a European option by 1-D quasi-Monte Carlo (Halton) integration.

    Converges to the exact Black-Scholes value much faster than pseudo-random
    Monte Carlo for the same ``n_points``. Deterministic (no seed needed).
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if n_points < 1:
        raise ValueError("n_points must be >= 1")
    if t == 0 or sigma == 0:
        return bsm_price(S, K, t, r, sigma, ot, b=b)

    drift = (b - 0.5 * sigma * sigma) * t
    vol = sigma * math.sqrt(t)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    total = 0.0
    for i in range(n_points):
        u = _radical_inverse(i + 1, 2)      # 1-D: base-2 van der Corput
        # Avoid the exact 0/1 endpoints where the inverse-CDF diverges.
        u = min(max(u, 1e-12), 1.0 - 1e-12)
        z = norm_ppf(u)
        sT = S * math.exp(drift + vol * z)
        total += max(sign * (sT - K), 0.0)
    return disc * total / n_points
