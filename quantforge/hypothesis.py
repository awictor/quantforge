"""Classical hypothesis tests built on the distribution CDFs.

Each test returns ``(statistic, p_value)`` computed against the relevant reference
distribution from :mod:`quantforge.distributions`:

- chi-square goodness-of-fit and test of independence (chi-square tail),
- one-way ANOVA (F tail),
- one-sample, paired and two-sample Student-t (pooled and Welch),
- the Mann-Whitney U rank-sum test (normal approximation), and
- the exact binomial test (binomial tail).

Pure standard library.
"""

import math

from .distributions import chi2_sf, f_cdf, binomial_pmf
from .student_t import t_cdf
from .mathfns import norm_cdf


def _mean(xs):
    return sum(xs) / len(xs)


def chi_square_gof_test(observed, expected=None):
    """Pearson chi-square goodness-of-fit test.

    Compares observed counts against ``expected`` (defaults to a uniform
    distribution over the categories). The statistic ``sum (O - E)^2 / E`` is
    referenced to a chi-square with ``k - 1`` degrees of freedom. Returns
    ``(statistic, p_value)``; a small p-value rejects the fit.
    """
    k = len(observed)
    if k < 2:
        raise ValueError("need at least two categories")
    if expected is None:
        total = sum(observed)
        expected = [total / k] * k
    if len(expected) != k:
        raise ValueError("observed and expected must have the same length")
    if any(e <= 0.0 for e in expected):
        raise ValueError("expected counts must be positive")
    stat = sum((o - e) ** 2 / e for o, e in zip(observed, expected))
    return stat, chi2_sf(stat, k - 1)


def chi_square_independence_test(table):
    """Chi-square test of independence for a contingency ``table`` (rows x cols).

    Uses the row/column marginals to form the expected counts under independence
    and references ``sum (O - E)^2 / E`` to a chi-square with
    ``(rows - 1)(cols - 1)`` degrees of freedom. Returns ``(statistic, p_value)``.
    """
    r = len(table)
    if r < 2 or any(len(row) != len(table[0]) for row in table):
        raise ValueError("table must be rectangular with at least two rows")
    c = len(table[0])
    if c < 2:
        raise ValueError("need at least two columns")
    row_sums = [sum(row) for row in table]
    col_sums = [sum(table[i][j] for i in range(r)) for j in range(c)]
    total = sum(row_sums)
    if total <= 0:
        raise ValueError("table total must be positive")
    stat = 0.0
    for i in range(r):
        for j in range(c):
            e = row_sums[i] * col_sums[j] / total
            if e <= 0.0:
                raise ValueError("zero expected count; collapse sparse categories")
            stat += (table[i][j] - e) ** 2 / e
    df = (r - 1) * (c - 1)
    return stat, chi2_sf(stat, df)


def one_way_anova(*groups):
    """One-way ANOVA F-test across two or more samples.

    Partitions the total variation into between-group and within-group sums of
    squares and forms ``F = MS_between / MS_within``, referenced to an F with
    ``(k - 1, N - k)`` degrees of freedom. Returns ``(F, p_value)``; a small
    p-value rejects equality of the group means.
    """
    groups = [list(g) for g in groups]
    k = len(groups)
    if k < 2:
        raise ValueError("need at least two groups")
    if any(len(g) < 1 for g in groups):
        raise ValueError("each group must be non-empty")
    n = sum(len(g) for g in groups)
    if n <= k:
        raise ValueError("need more observations than groups")
    grand = sum(sum(g) for g in groups) / n
    ss_between = sum(len(g) * (_mean(g) - grand) ** 2 for g in groups)
    ss_within = sum(sum((x - _mean(g)) ** 2 for x in g) for g in groups)
    df_between = k - 1
    df_within = n - k
    if ss_within == 0.0:
        raise ValueError("zero within-group variation")
    f = (ss_between / df_between) / (ss_within / df_within)
    return f, 1.0 - f_cdf(f, df_between, df_within)


def two_sample_t_test(a, b, equal_var=True):
    """Two-sample Student-t test of equal means, two-sided.

    With ``equal_var=True`` uses the pooled-variance t-test (``n_a + n_b - 2``
    degrees of freedom); with ``equal_var=False`` uses Welch's t-test with the
    Welch-Satterthwaite degrees of freedom. Returns ``(t, p_value)``.
    """
    a = list(a)
    b = list(b)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        raise ValueError("each sample needs at least two observations")
    ma, mb = _mean(a), _mean(b)
    va = sum((x - ma) ** 2 for x in a) / (na - 1)
    vb = sum((x - mb) ** 2 for x in b) / (nb - 1)
    if va == 0.0 and vb == 0.0:
        raise ValueError("both samples have zero variance")
    if equal_var:
        sp2 = ((na - 1) * va + (nb - 1) * vb) / (na + nb - 2)
        se = math.sqrt(sp2 * (1.0 / na + 1.0 / nb))
        df = na + nb - 2
    else:
        se = math.sqrt(va / na + vb / nb)
        df = (va / na + vb / nb) ** 2 / (
            (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    t = (ma - mb) / se
    # Two-sided p-value from the t CDF.
    p = 2.0 * (1.0 - t_cdf(abs(t), df))
    return t, p


def one_sample_t_test(sample, mu0=0.0):
    """One-sample two-sided Student-t test that the mean equals ``mu0``.

    ``t = (xbar - mu0) / (s / sqrt(n))`` on ``n - 1`` degrees of freedom. Returns
    ``(t, p_value)``.
    """
    x = list(sample)
    n = len(x)
    if n < 2:
        raise ValueError("need at least two observations")
    m = _mean(x)
    var = sum((v - m) ** 2 for v in x) / (n - 1)
    if var == 0.0:
        raise ValueError("zero-variance sample")
    t = (m - mu0) / math.sqrt(var / n)
    return t, 2.0 * (1.0 - t_cdf(abs(t), n - 1))


def paired_t_test(a, b):
    """Paired (dependent) two-sided Student-t test on the within-pair differences.

    Equivalent to a one-sample t-test of ``a[i] - b[i]`` against zero, on ``n - 1``
    degrees of freedom. Returns ``(t, p_value)``.
    """
    a = list(a)
    b = list(b)
    if len(a) != len(b):
        raise ValueError("paired samples must have the same length")
    if len(a) < 2:
        raise ValueError("need at least two pairs")
    diffs = [a[i] - b[i] for i in range(len(a))]
    return one_sample_t_test(diffs, 0.0)


def mann_whitney_u(a, b):
    """Mann-Whitney U rank-sum test (two-sided, normal approximation with ties).

    Ranks the pooled samples (average ranks for ties) and forms the smaller of the
    two U statistics; the p-value uses the normal approximation with a tie
    correction to the variance and a continuity correction. Returns
    ``(u, p_value)``, where ``u`` is ``min(U_a, U_b)``. A distribution-free
    alternative to the two-sample t when normality is doubtful.
    """
    a = list(a)
    b = list(b)
    na, nb = len(a), len(b)
    if na < 1 or nb < 1:
        raise ValueError("both samples must be non-empty")
    pooled = [(v, 0) for v in a] + [(v, 1) for v in b]
    pooled.sort(key=lambda t: t[0])
    n = na + nb
    # Average ranks for ties.
    ranks = [0.0] * n
    i = 0
    tie_term = 0.0
    while i < n:
        j = i
        while j + 1 < n and pooled[j + 1][0] == pooled[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0        # 1-based average rank
        for k in range(i, j + 1):
            ranks[k] = avg
        t = j - i + 1
        tie_term += t ** 3 - t
        i = j + 1
    rank_sum_a = sum(ranks[k] for k in range(n) if pooled[k][1] == 0)
    u_a = rank_sum_a - na * (na + 1) / 2.0
    u_b = na * nb - u_a
    u = min(u_a, u_b)
    mean_u = na * nb / 2.0
    var_u = na * nb / 12.0 * ((n + 1) - tie_term / (n * (n - 1)))
    if var_u <= 0.0:
        return u, 1.0
    z = (u - mean_u + 0.5) / math.sqrt(var_u)   # continuity correction toward mean
    p = 2.0 * norm_cdf(z)                        # u <= mean_u so z <= 0
    return u, min(1.0, p)


def binomial_test(k, n, prob=0.5, alternative="two-sided"):
    """Exact binomial test that the success probability equals ``prob``.

    ``alternative`` is ``"greater"`` (``P(X >= k)``), ``"less"`` (``P(X <= k)``), or
    ``"two-sided"`` (sum of all outcome probabilities no larger than the observed
    one). Returns ``(proportion, p_value)`` where ``proportion = k / n``.
    """
    if n < 1:
        raise ValueError("n must be at least 1")
    if not (0 <= k <= n):
        raise ValueError("require 0 <= k <= n")
    if not (0.0 <= prob <= 1.0):
        raise ValueError("prob must be in [0, 1]")
    pmf = [binomial_pmf(i, n, prob) for i in range(n + 1)]
    if alternative == "greater":
        p = sum(pmf[k:])
    elif alternative == "less":
        p = sum(pmf[:k + 1])
    elif alternative == "two-sided":
        threshold = pmf[k] * (1.0 + 1e-7)
        p = sum(pi for pi in pmf if pi <= threshold)
    else:
        raise ValueError("alternative must be 'greater', 'less' or 'two-sided'")
    return k / n, min(1.0, p)
