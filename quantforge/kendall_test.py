"""Tie-corrected Kendall's tau-b, Goodman-Kruskal gamma, and a significance test.

The tau-a in :mod:`quantforge.copula_stats` divides concordant-minus-discordant by
the raw pair count, so ties drag it below +/-1 even under a perfect monotone relation.
This module adds:

  * ``kendall_tau_b`` -- normalizes by the geometric mean of the untied-pair counts
    in each margin, so a perfect monotone fit reaches +/-1 despite ties;
  * ``goodman_kruskal_gamma`` -- ignores tied pairs entirely,
    ``(C - D) / (C + D)``;
  * ``kendall_tau_test`` -- the large-sample normal approximation for testing
    ``tau = 0`` (independence), returning tau-b, the z statistic and a two-sided
    p-value.

All pair counting is the direct ``O(n^2)`` sweep -- exact, no sort tricks. Pure
standard library.
"""

import math


def _counts(x, y):
    """Return (concordant, discordant, tied_x_only, tied_y_only, n)."""
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")
    c = d = tx = ty = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            if dx == 0 and dy == 0:
                continue                # tied on both -> counts in neither margin term
            if dx == 0:
                ty += 1
            elif dy == 0:
                tx += 1
            else:
                s = dx * dy
                if s > 0:
                    c += 1
                else:
                    d += 1
    return c, d, tx, ty, n


def kendall_tau_b(x, y):
    """Tie-corrected Kendall's tau-b.

    ``(C - D) / sqrt((C + D + Tx)(C + D + Ty))`` where ``C``/``D`` are concordant/
    discordant pairs and ``Tx``/``Ty`` are pairs tied only in ``x`` / only in ``y``.
    Reaches +/-1 for a perfect monotone relationship even when ties are present.
    """
    c, d, tx, ty, _ = _counts(x, y)
    denom = math.sqrt((c + d + tx) * (c + d + ty))
    if denom == 0:
        return 0.0
    return (c - d) / denom


def goodman_kruskal_gamma(x, y):
    """Goodman-Kruskal gamma: ``(C - D) / (C + D)``, ignoring all tied pairs."""
    c, d, _, _, _ = _counts(x, y)
    if c + d == 0:
        return 0.0
    return (c - d) / (c + d)


def kendall_tau_test(x, y):
    """Test ``tau = 0`` by the large-sample normal approximation.

    Returns a dict with ``tau_b``, the ``z`` statistic and the two-sided ``p_value``.
    The variance of the concordance statistic ``S = C - D`` under independence is
    ``n(n-1)(2n+5)/18`` (no tie correction to the variance -- adequate away from heavy
    ties); ``z = S / sqrt(Var S)``.
    """
    c, d, tx, ty, n = _counts(x, y)
    s = c - d
    var_s = n * (n - 1) * (2 * n + 5) / 18.0
    if var_s <= 0:
        z = 0.0
    else:
        z = s / math.sqrt(var_s)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return {"tau_b": kendall_tau_b(x, y), "z": z, "p_value": p}
