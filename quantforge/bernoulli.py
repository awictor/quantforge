"""Bernoulli numbers and Faulhaber's formula for sums of powers.

Bernoulli numbers ``B_n`` are the rational constants that appear in the closed form for
``sum_{k=1}^{m} k^p`` (Faulhaber's formula), the Euler-Maclaurin summation formula, and
the Taylor series of ``x / (e^x - 1)``. Computed exactly with :mod:`fractions`, so the
sum-of-powers polynomial is exact for any degree. Pure standard library.
"""

from fractions import Fraction
from math import comb, factorial

_BERNOULLI_CACHE = {0: Fraction(1)}


def bernoulli_number(n):
    """The ``n``-th Bernoulli number ``B_n`` as an exact :class:`fractions.Fraction`.

    Uses the convention ``B_1 = +1/2``. Computed by the recurrence
    ``sum_{k=0}^{n} C(n+1, k) B_k = 0`` and cached. ``B_n = 0`` for odd ``n > 1``.
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    if n in _BERNOULLI_CACHE:
        return _BERNOULLI_CACHE[n]
    for m in range(len(_BERNOULLI_CACHE), n + 1):
        if m == 1:
            _BERNOULLI_CACHE[1] = Fraction(1, 2)      # B_1 = +1/2 convention
            continue
        s = Fraction(0)
        for k in range(m):
            # The recurrence uses the B_1 = -1/2 convention; substitute it for k == 1.
            bk = Fraction(-1, 2) if k == 1 else _BERNOULLI_CACHE[k]
            s += comb(m + 1, k) * bk
        _BERNOULLI_CACHE[m] = -s / (m + 1)
    return _BERNOULLI_CACHE[n]


def faulhaber(m, p):
    """Sum ``1^p + 2^p + ... + m^p`` in closed form (Faulhaber), exact.

    Evaluates ``(1/(p+1)) sum_{j=0}^{p} C(p+1, j) B_j m^{p+1-j}`` with ``B_1 = +1/2``.
    Returns an ``int`` (the sum is always an integer). Handles ``p = 0`` (returns ``m``).
    """
    if m < 0:
        raise ValueError("m must be >= 0")
    if p < 0:
        raise ValueError("p must be >= 0")
    total = Fraction(0)
    for j in range(p + 1):
        total += comb(p + 1, j) * bernoulli_number(j) * Fraction(m) ** (p + 1 - j)
    total /= (p + 1)
    return int(total)


def bernoulli_sequence(n):
    """The Bernoulli numbers ``B_0 .. B_n`` as a list of :class:`Fraction`."""
    if n < 0:
        raise ValueError("n must be >= 0")
    return [bernoulli_number(k) for k in range(n + 1)]
