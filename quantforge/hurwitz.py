"""Hurwitz zeta and polygamma functions (Euler-Maclaurin summation).

The Hurwitz zeta ``hurwitz_zeta(s, a) = sum_{n>=0} (n + a)^-s`` (for ``s > 1``, ``a > 0``)
generalizes the Riemann zeta (``a = 1``) and is the natural home for the polygamma functions,
lattice sums, and Lerch-type series. As with :func:`quantforge.zeta.riemann_zeta`, the slow
tail is handled by Euler-Maclaurin with Bernoulli-number corrections. The polygamma
``polygamma(m, x) = d^m/dx^m psi(x)`` follows from the identity
``psi^(m)(x) = (-1)^(m+1) m! * hurwitz_zeta(m+1, x)`` for ``m >= 1``. Pure standard library.
"""

import math

from .bernoulli import bernoulli_number
from .special import digamma


def hurwitz_zeta(s, a, terms=20, corrections=10):
    """Hurwitz zeta ``zeta(s, a) = sum_{n>=0} (n + a)^-s`` for real ``s > 1``, ``a > 0``.

    Euler-Maclaurin summation: the first ``terms`` shifted terms are summed exactly, the tail
    is replaced by its integral, and ``corrections`` Bernoulli terms recover the remainder.
    With ``a = 1`` this reduces to the Riemann zeta.
    """
    if s <= 1.0:
        raise ValueError("hurwitz_zeta requires s > 1")
    if a <= 0.0:
        raise ValueError("hurwitz_zeta requires a > 0")
    N = terms
    # sum_{n=0}^{N-1} (n+a)^-s
    total = sum((n + a) ** (-s) for n in range(N))
    x = N + a
    total += x ** (1.0 - s) / (s - 1.0)   # integral tail
    total += 0.5 * x ** (-s)              # endpoint correction
    # Bernoulli corrections: sum_k B_2k/(2k)! * (s)_{2k-1} * x^{-s-2k+1}
    for k in range(1, corrections + 1):
        b2k = float(bernoulli_number(2 * k))
        prod = 1.0
        for j in range(2 * k - 1):
            prod *= (s + j)
        total += b2k / math.factorial(2 * k) * prod * x ** (-s - 2 * k + 1)
    return total


def polygamma(m, x):
    """Polygamma ``psi^(m)(x)``, the ``m``-th derivative of the digamma, for ``x > 0``.

    ``m = 0`` returns the digamma :func:`quantforge.special.digamma`. For ``m >= 1`` uses
    ``psi^(m)(x) = (-1)^(m+1) m! * hurwitz_zeta(m+1, x)``. ``polygamma(1, .)`` is the trigamma.
    """
    if m < 0:
        raise ValueError("polygamma order m must be >= 0")
    if x <= 0.0:
        raise ValueError("polygamma requires x > 0")
    if m == 0:
        return digamma(x)
    sign = 1.0 if (m + 1) % 2 == 0 else -1.0
    return sign * math.factorial(m) * hurwitz_zeta(m + 1, x)
