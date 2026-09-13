"""Dunn's test: pairwise rank comparisons after a Kruskal-Wallis rejection.

A significant Kruskal-Wallis test says *some* group differs but not *which*. Dunn's
(1964) test compares every pair of groups using the mean ranks from the single pooled
ranking (not re-ranked per pair), with the tie-corrected standard error

    SE_ij = sqrt( [ N(N+1)/12 - sum(t^3 - t)/(12(N-1)) ] * (1/n_i + 1/n_j) ),

and a z statistic ``(Rbar_i - Rbar_j) / SE_ij``. The raw two-sided p-values are then
adjusted for the ``k(k-1)/2`` comparisons (Bonferroni or Holm). Pure standard
library; builds on :mod:`quantforge.multiple_testing`.
"""

import math

from .multiple_testing import bonferroni, holm


def _midranks(values):
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def dunn_test(groups, adjust="holm"):
    """Dunn's post-hoc pairwise test after Kruskal-Wallis.

    ``groups`` is a sequence of samples. Returns a list of pairwise result dicts, one
    per unordered pair ``(i, j)``, each with ``groups`` ``(i, j)``, the ``z`` statistic,
    the raw ``p_value`` and the ``p_adjusted`` value. ``adjust`` is ``"holm"``
    (default), ``"bonferroni"`` or ``None`` for no correction. Uses one pooled ranking
    with the tie correction, so it is consistent with the Kruskal-Wallis H.
    """
    gs = [list(g) for g in groups if len(g) > 0]
    k = len(gs)
    if k < 2:
        raise ValueError("need at least 2 non-empty groups")

    pooled = [v for g in gs for v in g]
    n = len(pooled)
    ranks = _midranks(pooled)

    # Mean rank per group.
    mean_rank = []
    idx = 0
    for g in gs:
        s = sum(ranks[idx + i] for i in range(len(g)))
        idx += len(g)
        mean_rank.append(s / len(g))

    # Tie correction term sum(t^3 - t).
    tie_sum = 0.0
    svals = sorted(pooled)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and svals[j + 1] == svals[i]:
            j += 1
        t = j - i + 1
        if t > 1:
            tie_sum += t ** 3 - t
        i = j + 1

    sigma2_base = n * (n + 1) / 12.0 - tie_sum / (12.0 * (n - 1))

    pairs = []
    raw_p = []
    for i in range(k):
        for j in range(i + 1, k):
            se = math.sqrt(sigma2_base * (1.0 / len(gs[i]) + 1.0 / len(gs[j])))
            z = (mean_rank[i] - mean_rank[j]) / se if se > 0 else 0.0
            p = math.erfc(abs(z) / math.sqrt(2.0))
            pairs.append((i, j, z, p))
            raw_p.append(p)

    if adjust == "bonferroni":
        adj = bonferroni(raw_p)
    elif adjust == "holm":
        adj = holm(raw_p)
    elif adjust is None:
        adj = raw_p
    else:
        raise ValueError("adjust must be 'holm', 'bonferroni' or None")

    return [
        {"groups": (i, j), "z": z, "p_value": p, "p_adjusted": adj[idx]}
        for idx, (i, j, z, p) in enumerate(pairs)
    ]
