"""Tail and shape metrics of the smile-implied risk-neutral density.

Given the Breeden-Litzenberger density ``g`` (see
:mod:`quantforge.rnd`), these are model-free summaries the desk reads straight
off the smile:

  * ``tail_probability`` -- risk-neutral ``Q(S_T < L)`` or ``Q(S_T > U)`` (the
    price of a cash-or-nothing digital, undiscounted);
  * ``density_entropy`` -- the differential entropy of the terminal spot's
    density, a spread/uncertainty measure;
  * ``expected_shortfall`` -- the risk-neutral conditional expectation of the
    spot in a tail, ``E^Q[S_T | S_T < L]`` (or above ``U``).

All integrate the density grid by the trapezoidal rule. Pure standard library.
"""

import math

from .rnd import density_grid_from_smile


def _grid(S0, t, r, vol_fn, q, n, width):
    ks, dens = density_grid_from_smile(S0, t, r, vol_fn, q=q, n=n, width=width)
    dens = [max(d, 0.0) for d in dens]           # clip tiny negative wiggles
    # Renormalise to unit mass (finite grid truncates a little).
    dK = ks[1] - ks[0]
    mass = sum(dens) * dK
    if mass > 0:
        dens = [d / mass for d in dens]
    return ks, dens


def tail_probability(S0, t, r, vol_fn, level, lower=True, q=0.0, n=600,
                     width=8.0):
    """Risk-neutral tail probability ``Q(S_T < level)`` (or ``> level``).

    Equals the undiscounted price of a cash-or-nothing binary struck at
    ``level``. ``lower=True`` returns the downside probability.
    """
    ks, dens = _grid(S0, t, r, vol_fn, q, n, width)
    total = 0.0
    for i in range(len(ks) - 1):
        k0, k1 = ks[i], ks[i + 1]
        inside0 = (k0 < level) if lower else (k0 > level)
        inside1 = (k1 < level) if lower else (k1 > level)
        if inside0 and inside1:
            total += 0.5 * (dens[i] + dens[i + 1]) * (k1 - k0)
    return total


def density_entropy(S0, t, r, vol_fn, q=0.0, n=600, width=8.0):
    """Differential entropy ``-integral g ln g dK`` of the terminal-spot density."""
    ks, dens = _grid(S0, t, r, vol_fn, q, n, width)
    total = 0.0
    for i in range(len(ks) - 1):
        k0, k1 = ks[i], ks[i + 1]
        h0 = -dens[i] * math.log(dens[i]) if dens[i] > 1e-300 else 0.0
        h1 = -dens[i + 1] * math.log(dens[i + 1]) if dens[i + 1] > 1e-300 else 0.0
        total += 0.5 * (h0 + h1) * (k1 - k0)
    return total


def expected_shortfall(S0, t, r, vol_fn, level, lower=True, q=0.0, n=600,
                       width=8.0):
    """Risk-neutral tail mean ``E^Q[S_T | S_T < level]`` (or ``> level``).

    Returns the conditional expectation of the terminal spot in the tail beyond
    ``level``; ``nan`` if that tail has zero probability.
    """
    ks, dens = _grid(S0, t, r, vol_fn, q, n, width)
    num = 0.0
    den = 0.0
    for i in range(len(ks) - 1):
        k0, k1 = ks[i], ks[i + 1]
        inside0 = (k0 < level) if lower else (k0 > level)
        inside1 = (k1 < level) if lower else (k1 > level)
        if inside0 and inside1:
            p = 0.5 * (dens[i] + dens[i + 1]) * (k1 - k0)
            mid = 0.5 * (k0 + k1)
            den += p
            num += mid * p
    return num / den if den > 1e-14 else float("nan")
