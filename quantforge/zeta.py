"""Riemann zeta and Dirichlet eta functions (Euler-Maclaurin summation).

The Riemann zeta ``zeta(s) = sum_{n>=1} n^-s`` (for ``s > 1``) is the backbone of analytic
number theory and appears in Casimir energy, prime counting, and regularization. Direct
summation converges slowly, so this uses Euler-Maclaurin: sum the first ``N`` terms exactly,
approximate the tail by an integral, and correct with Bernoulli-number terms -- accurate for
any real ``s > 1`` (and, via the eta relation, for ``0 < s < 1``). The Dirichlet eta
``eta(s) = sum (-1)^{n-1} n^-s = (1 - 2^{1-s}) zeta(s)`` converges for ``s > 0``. Pure
standard library.
"""

import math

from .bernoulli import bernoulli_number


def riemann_zeta(s, terms=20, corrections=10):
    """Riemann zeta ``zeta(s)`` for real ``s > 0``, ``s != 1`` (Euler-Maclaurin).

    For ``s > 1`` uses Euler-Maclaurin directly; for ``0 < s < 1`` uses the Dirichlet-eta
    relation ``zeta(s) = eta(s) / (1 - 2^{1-s})``. Raises at the pole ``s = 1``.
    """
    if abs(s - 1.0) < 1e-12:
        raise ValueError("zeta has a pole at s = 1")
    if s < 1.0:
        # analytic continuation on (0, 1) via eta (which converges there)
        return dirichlet_eta(s) / (1.0 - 2.0 ** (1.0 - s))
    return _zeta_em(s, terms, corrections)


def _zeta_em(s, N, K):
    # Euler-Maclaurin: sum_{n=1}^{N-1} n^-s + N^-s/2 + N^{1-s}/(s-1)
    #   + sum_{k=1}^{K} B_{2k}/(2k)! * (falling factorial) * N^{-s-2k+1}
    total = sum(n ** (-s) for n in range(1, N))
    total += 0.5 * N ** (-s)
    total += N ** (1 - s) / (s - 1)
    # correction terms
    term_coeff = s          # (s)(s+1)...(s+2k-2) built up
    for k in range(1, K + 1):
        b2k = float(bernoulli_number(2 * k))
        # product s(s+1)...(s+2k-2) has 2k-1 factors
        prod = 1.0
        for j in range(2 * k - 1):
            prod *= (s + j)
        total += b2k / math.factorial(2 * k) * prod * N ** (-s - 2 * k + 1)
    return total


def dirichlet_eta(s):
    """Dirichlet eta ``eta(s) = sum (-1)^{n-1} n^-s = (1 - 2^{1-s}) zeta(s)`` for ``s > 0``.

    Converges (conditionally) for all ``s > 0``; accelerated here by van Wijngaarden /
    alternating-series transformation for robustness near ``s -> 0``.
    """
    if s <= 0:
        raise ValueError("dirichlet_eta requires s > 0")
    if s > 1.0:
        return (1.0 - 2.0 ** (1.0 - s)) * _zeta_em(s, 20, 10)
    # alternating-series acceleration (Cohen-Villegas-Zagier) for 0 < s <= 1
    n = 40
    d = (3 + math.sqrt(8)) ** n
    d = (d + 1 / d) / 2
    b = -1.0
    c = -d
    total = 0.0
    for k in range(n):
        c = b - c
        total += c * ((k + 1) ** (-s))
        b = (k + n) * (k - n) * b / ((k + 0.5) * (k + 1))
    return total / d
