"""Anderson-Darling k-sample test for a common distribution.

Tests the null that ``k`` independent samples are drawn from the *same* (unspecified)
distribution -- the nonparametric analogue of one-way ANOVA, but sensitive to
differences anywhere in the distribution (location, scale, shape), not just the mean.
It compares each sample's empirical CDF to the pooled CDF via the Anderson-Darling
weighting (Scholz-Stephens 1987), which emphasizes the tails. Returns a standardized
statistic and an approximate p-value from the tabulated null. Pure standard library.
"""

import math


def anderson_darling_ksample(*samples):
    """Scholz-Stephens k-sample Anderson-Darling test.

    Pass two or more samples as separate sequence arguments. Returns a dict with the
    raw statistic ``a2k``, the ``standardized`` statistic ``(A2k - (k-1)) / sqrt(var)``,
    and an approximate ``p_value``. A small p-value rejects the null that all samples
    share one distribution.
    """
    samples = [sorted(s) for s in samples if len(s) > 0]
    k = len(samples)
    if k < 2:
        raise ValueError("need at least 2 non-empty samples")
    ns = [len(s) for s in samples]
    N = sum(ns)
    pooled = sorted(v for s in samples for v in s)

    # Distinct pooled values and their multiplicities (mid-rank handling of ties).
    distinct = []
    mult = []
    for v in pooled:
        if distinct and distinct[-1] == v:
            mult[-1] += 1
        else:
            distinct.append(v)
            mult.append(1)
    L = len(distinct)
    if L < 2:
        raise ValueError("all values identical; test undefined")

    # Scholz-Stephens A2akN statistic (version with ties, their eq. 7-8).
    # B_j = cumulative pooled count up to and including value j (using mid-ranks).
    a2k = 0.0
    for i in range(k):
        s = samples[i]
        ni = ns[i]
        # Counts of sample i at each distinct value.
        cnt_i = {}
        for v in s:
            cnt_i[v] = cnt_i.get(v, 0) + 1
        cum_pool = 0      # l_j cumulative pooled before j
        cum_i = 0         # cumulative in sample i before j
        inner = 0.0
        for j in range(L):
            lj = mult[j]
            fij = cnt_i.get(distinct[j], 0)
            # Mid-point cumulative counts (Bj uses count strictly before + half of tie).
            Bj = cum_pool + lj / 2.0
            Mij = cum_i + fij / 2.0
            denom = Bj * (N - Bj) - N * lj / 4.0
            if denom > 0:
                inner += lj * (N * Mij - ni * Bj) ** 2 / denom
            cum_pool += lj
            cum_i += fij
        a2k += inner / ni
    a2k *= (N - 1) / (N * N)

    # Null mean is k-1; variance from Scholz-Stephens eq. 4.
    H = sum(1.0 / n for n in ns)
    h = sum(1.0 / i for i in range(1, N))
    g = 0.0
    for i in range(1, N - 1):
        for j in range(i + 1, N):
            g += 1.0 / ((N - i) * j)
    a = (4 * g - 6) * (k - 1) + (10 - 6 * g) * H
    b = (2 * g - 4) * k * k + 8 * h * k + (2 * g - 14 * h - 4) * H - 8 * h + 4 * g - 6
    c = (6 * h + 2 * g - 2) * k * k + (4 * h - 4 * g + 6) * k + (2 * h - 6) * H + 4 * h
    d = (2 * h + 6) * k * k - 4 * h * k
    var = (a * N ** 3 + b * N ** 2 + c * N + d) / ((N - 1) * (N - 2) * (N - 3))
    std = (a2k - (k - 1)) / math.sqrt(var) if var > 0 else 0.0

    # Approximate p-value from the Scholz-Stephens tabulated upper tail (interpolated).
    p = _ad_ksample_pvalue(std)
    return {"a2k": a2k, "standardized": std, "p_value": p}


def _ad_ksample_pvalue(tm):
    """Interpolate the Scholz-Stephens p-value from the standardized statistic."""
    # Tabulated critical values (m = k-1 large-sample limit), Scholz-Stephens table 1.
    t = [0.325, 1.226, 1.960, 2.719, 3.752, 4.592, 6.546]      # standardized crit
    p = [0.25, 0.10, 0.05, 0.025, 0.01, 0.005, 0.001]
    if tm <= t[0]:
        return 1.0 if tm < t[0] else 0.25
    if tm >= t[-1]:
        return 0.001
    for i in range(len(t) - 1):
        if t[i] <= tm <= t[i + 1]:
            # Log-linear interpolation in p.
            frac = (tm - t[i]) / (t[i + 1] - t[i])
            return math.exp(math.log(p[i]) + frac * (math.log(p[i + 1]) - math.log(p[i])))
    return 0.001
