"""Rank-based dependence: Kendall's tau, Spearman's rho, empirical copula.

Linear (Pearson) correlation only captures linear co-movement and is not invariant
to monotone transforms. Rank-based measures fix both: they depend only on the
*ordering* of the data, so they are invariant to any strictly increasing marginal
transform and capture monotone (not just linear) dependence -- the natural
language of copulas.

  * ``kendall_tau``   -- (concordant - discordant) pairs / total pairs, in [-1, 1].
  * ``spearman_rho``  -- Pearson correlation of the ranks, in [-1, 1].
  * ``pseudo_observations`` -- map each margin to its scaled rank in (0, 1); the
    joint of these is the empirical copula, stripped of the marginals.

Pure standard library.
"""


def _ranks(x):
    """Average ranks (1-based), ties share the mean of their positions."""
    n = len(x)
    order = sorted(range(n), key=lambda i: x[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and x[order[j + 1]] == x[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0        # average of positions i..j (1-based)
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def kendall_tau(x, y):
    """Kendall's rank correlation tau-a between paired samples.

    Counts concordant minus discordant pairs over all ``n(n-1)/2`` pairs. Equals
    +1 for a strictly increasing relationship, -1 for strictly decreasing, and ~0
    under independence. This is the tau-a variant (no tie correction).
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have the same length")
    if n < 2:
        raise ValueError("need at least 2 points")
    concordant = 0
    discordant = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            s = dx * dy
            if s > 0.0:
                concordant += 1
            elif s < 0.0:
                discordant += 1
    total = n * (n - 1) // 2
    return (concordant - discordant) / total


def spearman_rho(x, y):
    """Spearman's rank correlation: the Pearson correlation of the ranks.

    Invariant to any monotone transform of either margin; +1/-1 for a perfectly
    monotone relationship, ~0 under independence.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have the same length")
    if n < 2:
        raise ValueError("need at least 2 points")
    rx = _ranks(x)
    ry = _ranks(y)
    mx = sum(rx) / n
    my = sum(ry) / n
    sxy = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    sxx = sum((rx[i] - mx) ** 2 for i in range(n))
    syy = sum((ry[i] - my) ** 2 for i in range(n))
    if sxx <= 0.0 or syy <= 0.0:
        raise ValueError("a variable has no rank variation (all ties)")
    return sxy / (sxx * syy) ** 0.5


def pseudo_observations(x):
    """Empirical-copula pseudo-observations: scaled ranks in (0, 1).

    Maps each value to ``rank / (n + 1)``, the standard normalization that keeps
    the transformed sample strictly inside the open unit interval. Applying this
    to each margin and viewing the joint gives the empirical copula, free of the
    marginal distributions.
    """
    n = len(x)
    if n < 1:
        raise ValueError("need at least 1 point")
    r = _ranks(x)
    return [ri / (n + 1.0) for ri in r]
