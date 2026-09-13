"""Jonckheere-Terpstra test for an ordered (monotone) trend across k groups.

Kruskal-Wallis asks whether the groups differ *at all*; when the groups have a
natural order (dose levels, time buckets, rating tiers) and you expect the response
to move monotonically, the Jonckheere-Terpstra test is far more powerful. For every
ordered pair of groups ``(i < j)`` it counts the Mann-Whitney concordances -- how
often an observation in the later group exceeds one in the earlier -- and sums them:

    J = sum_{i < j} #{ (x in group_i, y in group_j) : y > x } + 0.5 * ties.

Under the null of no trend ``J`` is approximately normal with a known mean and
variance (with a tie correction); a large positive ``z`` signals an increasing trend
across the group order, a large negative one a decreasing trend. Pure standard
library.
"""

import math


def _mann_whitney_count(a, b):
    """Sum over pairs of 1 if b_j > a_i, 0.5 if tied (the U-count of b over a)."""
    count = 0.0
    for x in a:
        for y in b:
            if y > x:
                count += 1.0
            elif y == x:
                count += 0.5
    return count


def jonckheere_terpstra_test(groups):
    """Jonckheere-Terpstra trend test across ordered ``groups``.

    ``groups`` is a sequence of samples given in the hypothesized order (e.g. from
    lowest dose to highest). Returns a dict with the ``statistic`` J, its null ``mean``
    and ``variance`` (tie-corrected), the ``z`` normal approximation and the two-sided
    ``p_value``. A positive ``z`` indicates an increasing trend across the group order.
    """
    gs = [list(g) for g in groups if len(g) > 0]
    k = len(gs)
    if k < 2:
        raise ValueError("need at least 2 non-empty groups")

    j_stat = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            j_stat += _mann_whitney_count(gs[i], gs[j])

    n = sum(len(g) for g in gs)
    sizes = [len(g) for g in gs]
    mean = (n * n - sum(s * s for s in sizes)) / 4.0

    # Tie-corrected variance (Lehmann). Pool all values for the tie term.
    pooled = [v for g in gs for v in g]
    tie_term_t = 0.0
    svals = sorted(pooled)
    i = 0
    while i < n:
        r = i
        while r + 1 < n and svals[r + 1] == svals[i]:
            r += 1
        t = r - i + 1
        if t > 1:
            tie_term_t += t * (t - 1) * (2 * t + 5)
        i = r + 1

    sum_n = n * (n - 1) * (2 * n + 5)
    sum_ni = sum(s * (s - 1) * (2 * s + 5) for s in sizes)

    var = (sum_n - sum_ni - tie_term_t) / 72.0
    # Second-order tie terms (small; included for correctness with heavy ties).
    if tie_term_t > 0 and n > 2:
        sum_ni2 = sum(s * (s - 1) * (s - 2) for s in sizes)
        tie2 = sum(_t3(svals))  # sum t(t-1)(t-2)
        var += sum_ni2 * tie2 / (36.0 * n * (n - 1) * (n - 2))
        sum_ni_pair = sum(s * (s - 1) for s in sizes)
        tie_pair = sum(_t2(svals))  # sum t(t-1)
        var += sum_ni_pair * tie_pair / (8.0 * n * (n - 1))

    if var <= 0:
        z = 0.0
    else:
        z = (j_stat - mean) / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return {"statistic": j_stat, "mean": mean, "variance": var, "z": z, "p_value": p}


def _t2(svals):
    """Yield t(t-1) for each tie block in a sorted list."""
    n = len(svals)
    i = 0
    while i < n:
        r = i
        while r + 1 < n and svals[r + 1] == svals[i]:
            r += 1
        t = r - i + 1
        yield t * (t - 1)
        i = r + 1


def _t3(svals):
    """Yield t(t-1)(t-2) for each tie block in a sorted list."""
    n = len(svals)
    i = 0
    while i < n:
        r = i
        while r + 1 < n and svals[r + 1] == svals[i]:
            r += 1
        t = r - i + 1
        yield t * (t - 1) * (t - 2)
        i = r + 1
