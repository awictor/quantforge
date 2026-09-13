"""Kruskal-Wallis and Friedman rank tests: nonparametric k-sample comparisons.

Two distribution-free extensions of the two-sample rank tests to more than two groups:

  * ``kruskal_wallis_test`` -- the k-sample analogue of the Mann-Whitney / one-way
    ANOVA. Pool all observations, rank them, and compare each group's mean rank to the
    overall mean rank; the statistic ``H`` is chi-square with ``k - 1`` degrees of
    freedom under equal distributions. On two groups it agrees with the Mann-Whitney
    normal approximation.
  * ``friedman_test`` -- the repeated-measures analogue (a nonparametric two-way ANOVA
    with one observation per cell). Rank *within* each block (row) across the ``k``
    treatments, then compare the treatments' rank sums; the statistic is chi-square
    with ``k - 1`` degrees of freedom.

Both apply the standard tie correction and return a p-value from the chi-square upper
tail. Pure standard library.
"""

from .distributions import chi2_cdf


def _midranks(values):
    """1-based midranks of a flat list (ties averaged)."""
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


def kruskal_wallis_test(*groups):
    """Kruskal-Wallis H test that ``k`` groups share a distribution.

    Pass each group as a separate sequence argument. Returns a dict with the
    tie-corrected ``statistic`` H, the ``df`` (``k - 1``) and the chi-square upper-tail
    ``p_value``. A small p-value rejects the null of equal distributions (specifically,
    equal medians for similarly-shaped groups).
    """
    groups = [list(g) for g in groups if len(g) > 0]
    k = len(groups)
    if k < 2:
        raise ValueError("need at least 2 non-empty groups")
    pooled = [v for g in groups for v in g]
    n = len(pooled)
    ranks = _midranks(pooled)

    # Rank sums per group, walking the pooled rank list in group order.
    idx = 0
    h_sum = 0.0
    for g in groups:
        rsum = sum(ranks[idx + i] for i in range(len(g)))
        idx += len(g)
        h_sum += rsum * rsum / len(g)

    h = 12.0 / (n * (n + 1)) * h_sum - 3.0 * (n + 1)

    # Tie correction: divide by 1 - sum(t^3 - t) / (N^3 - N).
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
    if tie_sum > 0 and n ** 3 - n > 0:
        h /= 1.0 - tie_sum / (n ** 3 - n)

    df = k - 1
    p = 1.0 - chi2_cdf(h, df)
    return {"statistic": h, "df": df, "p_value": p}


def friedman_test(blocks):
    """Friedman test for ``k`` related treatments over ``b`` blocks.

    ``blocks`` is a sequence of rows, each a length-``k`` sequence giving one block's
    measurements across the treatments (e.g. one subject rated under every condition).
    Ranks within each block, then compares treatment rank sums. Returns a dict with
    the ``statistic``, ``df`` (``k - 1``) and the chi-square ``p_value``.
    """
    rows = [list(r) for r in blocks]
    b = len(rows)
    if b < 2:
        raise ValueError("need at least 2 blocks")
    k = len(rows[0])
    if k < 2:
        raise ValueError("need at least 2 treatments")
    if any(len(r) != k for r in rows):
        raise ValueError("all blocks must have the same number of treatments")

    col_rank_sum = [0.0] * k
    for r in rows:
        rr = _midranks(r)
        for c in range(k):
            col_rank_sum[c] += rr[c]

    stat = 12.0 / (b * k * (k + 1)) * sum(rs * rs for rs in col_rank_sum) \
        - 3.0 * b * (k + 1)
    df = k - 1
    p = 1.0 - chi2_cdf(stat, df)
    return {"statistic": stat, "df": df, "p_value": p}
