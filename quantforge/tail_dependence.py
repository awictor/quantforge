"""Empirical tail dependence and exceedance correlation.

Correlation measures average co-movement; it says nothing about whether two assets
crash *together*. Tail dependence does: the empirical upper tail-dependence
coefficient estimates ``P(U > q | V > q)`` as the threshold ``q -> 1`` (both
margins on their rank/uniform scale), and the lower coefficient the analogous
``P(U <= q | V <= q)`` as ``q -> 0``. A value near 1 means extremes cluster; near 0
means joint extremes are no more likely than under independence.

Exceedance correlation is the correlation computed only on observations where both
margins breach a quantile -- the "correlations rise in a crash" effect. These are
estimated directly from ranks (empirical copula), free of the marginals. Pure
standard library.
"""


def _ranks_uniform(x):
    """Scaled ranks in (0, 1): rank / (n + 1). Ties share their average rank."""
    n = len(x)
    order = sorted(range(n), key=lambda i: x[i])
    u = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and x[order[j + 1]] == x[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            u[order[k]] = avg / (n + 1.0)
        i = j + 1
    return u


def upper_tail_dependence(x, y, q=0.95):
    """Empirical upper tail-dependence coefficient at threshold ``q``.

    Estimates ``P(U > q | V > q)`` on the rank/uniform scale: of the points whose
    ``y`` rank exceeds ``q``, the fraction whose ``x`` rank also exceeds ``q``.
    Near 1 means the two crash/spike together; near ``1 - q``-scaled independence
    means they do not. Returns 0 when no point exceeds the threshold in ``y``.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if not (0.0 < q < 1.0):
        raise ValueError("q must be in (0, 1)")
    u = _ranks_uniform(x)
    v = _ranks_uniform(y)
    cond = [i for i in range(len(u)) if v[i] > q]
    if not cond:
        return 0.0
    joint = sum(1 for i in cond if u[i] > q)
    return joint / len(cond)


def lower_tail_dependence(x, y, q=0.05):
    """Empirical lower tail-dependence coefficient at threshold ``q``.

    Estimates ``P(U <= q | V <= q)`` on the rank/uniform scale. Near 1 means joint
    downside extremes cluster. Returns 0 when no point falls below the threshold
    in ``y``.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if not (0.0 < q < 1.0):
        raise ValueError("q must be in (0, 1)")
    u = _ranks_uniform(x)
    v = _ranks_uniform(y)
    cond = [i for i in range(len(u)) if v[i] <= q]
    if not cond:
        return 0.0
    joint = sum(1 for i in cond if u[i] <= q)
    return joint / len(cond)


def exceedance_correlation(x, y, q=0.9, tail="upper"):
    """Pearson correlation computed only on joint-tail observations.

    Selects the points where both margins breach the ``q`` quantile (upper tail)
    or fall below the ``1 - q`` quantile (lower tail), then returns the ordinary
    correlation of ``x`` and ``y`` on that subset -- the "correlations rise in the
    tails" diagnostic. Raises if fewer than two joint-tail points exist.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if not (0.0 < q < 1.0):
        raise ValueError("q must be in (0, 1)")
    u = _ranks_uniform(x)
    v = _ranks_uniform(y)
    if tail == "upper":
        idx = [i for i in range(len(u)) if u[i] > q and v[i] > q]
    elif tail == "lower":
        lo = 1.0 - q
        idx = [i for i in range(len(u)) if u[i] <= lo and v[i] <= lo]
    else:
        raise ValueError("tail must be 'upper' or 'lower'")
    if len(idx) < 2:
        raise ValueError("not enough joint-tail observations")
    xs = [x[i] for i in idx]
    ys = [y[i] for i in idx]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    sxy = sum((xs[k] - mx) * (ys[k] - my) for k in range(len(xs)))
    sxx = sum((xs[k] - mx) ** 2 for k in range(len(xs)))
    syy = sum((ys[k] - my) ** 2 for k in range(len(ys)))
    if sxx <= 0.0 or syy <= 0.0:
        raise ValueError("no variation in a tail subset")
    return sxy / (sxx * syy) ** 0.5
