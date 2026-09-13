"""Fisher's exact test for a 2x2 contingency table.

For a 2x2 table with small counts the chi-square approximation is unreliable; Fisher's
exact test computes the p-value exactly from the hypergeometric distribution.
Conditioning on the row and column margins, the probability of a table with top-left
cell ``a`` is

    P(a) = C(r1, a) C(r2, c1 - a) / C(n, c1),

where ``r1, r2`` are the row sums, ``c1`` the first column sum and ``n`` the total. The
two-sided p-value sums the probabilities of all tables no more likely than the
observed one. Also returns the sample odds ratio. Pure standard library.
"""

from math import comb


def _hypergeom_pmf(a, r1, r2, c1):
    n = r1 + r2
    return comb(r1, a) * comb(r2, c1 - a) / comb(n, c1)


def fisher_exact_test(table, alternative="two-sided"):
    """Fisher's exact test for a 2x2 ``table`` ``[[a, b], [c, d]]``.

    ``alternative`` is ``"two-sided"`` (default), ``"greater"`` or ``"less"`` (one-sided
    on the odds ratio). Returns a dict with the sample ``odds_ratio``
    (``a d / (b c)``, ``inf`` if a denominator is zero) and the ``p_value``. The
    two-sided p-value sums the probabilities of every table (given the margins) no more
    probable than the observed one.
    """
    if len(table) != 2 or any(len(row) != 2 for row in table):
        raise ValueError("table must be 2x2")
    a, b = table[0]
    c, d = table[1]
    if min(a, b, c, d) < 0:
        raise ValueError("counts must be non-negative")

    r1 = a + b
    r2 = c + d
    c1 = a + c
    n = r1 + r2
    if n == 0:
        raise ValueError("table is empty")

    # Odds ratio.
    if b * c == 0:
        odds = float("inf") if a * d > 0 else float("nan")
    else:
        odds = (a * d) / (b * c)

    # Support of the top-left cell given the margins.
    lo = max(0, c1 - r2)
    hi = min(r1, c1)

    p_obs = _hypergeom_pmf(a, r1, r2, c1)

    if alternative == "two-sided":
        tol = p_obs * (1.0 + 1e-9)
        p = sum(_hypergeom_pmf(k, r1, r2, c1)
                for k in range(lo, hi + 1)
                if _hypergeom_pmf(k, r1, r2, c1) <= tol)
    elif alternative == "greater":
        # Larger a -> larger odds ratio.
        p = sum(_hypergeom_pmf(k, r1, r2, c1) for k in range(a, hi + 1))
    elif alternative == "less":
        p = sum(_hypergeom_pmf(k, r1, r2, c1) for k in range(lo, a + 1))
    else:
        raise ValueError("alternative must be 'two-sided', 'greater' or 'less'")

    return {"odds_ratio": odds, "p_value": min(p, 1.0)}
