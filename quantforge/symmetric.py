"""Symmetric functions and Newton's identities.

Three ways to describe the same information about a multiset of values ``x_1..x_n``:

  * *elementary symmetric* polynomials ``e_k`` = the sum of all products of ``k`` distinct
    values (``e_0 = 1``), which are the signed coefficients of ``prod (t - x_i)``;
  * *power sums* ``p_k = sum x_i^k``;
  * the monic polynomial coefficients themselves (Vieta's formulas).

Newton's identities convert between the ``e_k`` and the ``p_k``. This module computes each
directly and converts among them, exactly for integer or rational (`fractions.Fraction`)
input. Pure standard library.
"""

from fractions import Fraction


def elementary_symmetric(values):
    """Return ``[e_0, e_1, ..., e_n]`` for ``values`` (``e_0 = 1``).

    ``e_k`` is the sum of all products of ``k`` distinct entries. Built by the standard
    ``O(n^2)`` DP that multiplies in one value at a time.
    """
    e = [1]
    for x in values:
        # multiply the generating polynomial by (1 + x t): e_k += x * e_{k-1}
        new = e + [0]
        for k in range(len(e), 0, -1):
            new[k] = new[k] + x * e[k - 1]
        e = new
    return e


def power_sums(values, kmax):
    """Return ``[p_0, p_1, ..., p_kmax]`` where ``p_k = sum x_i^k`` (``p_0 = n``)."""
    if kmax < 0:
        raise ValueError("kmax must be non-negative")
    n = len(values)
    out = [n]
    for k in range(1, kmax + 1):
        out.append(sum(x ** k for x in values))
    return out


def power_to_elementary(p):
    """Convert power sums ``[p_1, ..., p_n]`` to ``[e_0, ..., e_n]`` (Newton's identities).

    ``k e_k = sum_{i=1}^{k} (-1)^(i-1) e_{k-i} p_i``. Uses exact `Fraction` arithmetic so the
    division by ``k`` is exact for integer power sums.
    """
    n = len(p)
    e = [Fraction(1)]
    for k in range(1, n + 1):
        acc = Fraction(0)
        for i in range(1, k + 1):
            acc += (-1) ** (i - 1) * e[k - i] * p[i - 1]
        e.append(acc / k)
    return e


def elementary_to_power(e, kmax=None):
    """Convert ``[e_0, ..., e_n]`` to power sums ``[p_1, ..., p_kmax]`` (Newton's identities).

    ``p_k = sum_{i=1}^{k-1} (-1)^(i-1) e_i p_{k-i} + (-1)^(k-1) k e_k`` (with ``e_k = 0``
    for ``k > n``). Defaults ``kmax`` to ``n``.
    """
    n = len(e) - 1
    if kmax is None:
        kmax = n

    def ek(k):
        return e[k] if k <= n else 0

    p = []
    for k in range(1, kmax + 1):
        acc = 0
        for i in range(1, k):
            acc += (-1) ** (i - 1) * ek(i) * p[k - i - 1]
        acc += (-1) ** (k - 1) * k * ek(k)
        p.append(acc)
    return p


def poly_from_roots(roots):
    """Monic polynomial coefficients (highest-degree first) of ``prod (t - r)``.

    ``coeffs[k] = (-1)^k e_k``. Matches the coefficient order of
    :func:`quantforge.polyroots.polynomial_roots`.
    """
    e = elementary_symmetric(roots)
    return [(-1) ** k * e[k] for k in range(len(e))]
