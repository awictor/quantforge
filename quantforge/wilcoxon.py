"""Wilcoxon signed-rank test and the sign test for paired / one-sample data.

Distribution-free alternatives to the paired t-test. Both test whether a sample (or
the differences of a matched pair) is centered at a hypothesized median:

  * ``wilcoxon_signed_rank_test`` -- ranks the absolute deviations, sums the ranks of
    the positive ones, and compares to its null distribution. Uses the signed values,
    so it is more powerful than the sign test when the distribution is roughly
    symmetric; it pairs with the Hodges-Lehmann location estimator.
  * ``sign_test`` -- counts how many observations exceed the hypothesized median and
    tests that count against Binomial(n, 1/2). It uses only the *signs*, so it makes
    no symmetry assumption and is maximally robust (but less powerful).

Both return a normal-approximation two-sided p-value (with tie/zero handling);
the sign test also exposes the exact binomial p-value. Pure standard library.
"""

import math

from .distributions import binomial_cdf


def _erfc(z):
    return math.erfc(z)


def wilcoxon_signed_rank_test(x, mu0=0.0, y=None):
    """Wilcoxon signed-rank test that the (paired) sample is centered at ``mu0``.

    If ``y`` is given the test runs on the paired differences ``x_i - y_i``; otherwise
    on ``x_i - mu0``. Zero differences are dropped (Wilcoxon's convention); tied
    absolute values receive average ranks. Returns a dict with ``statistic`` W (the
    positive-rank sum), the ``z`` normal approximation (continuity-corrected, with the
    tie correction to the variance) and the two-sided ``p_value``.
    """
    if y is not None:
        if len(x) != len(y):
            raise ValueError("x and y must have equal length")
        diffs = [x[i] - y[i] for i in range(len(x))]
    else:
        diffs = [xi - mu0 for xi in x]

    nz = [d for d in diffs if d != 0.0]
    n = len(nz)
    if n == 0:
        return {"statistic": 0.0, "z": 0.0, "p_value": 1.0, "n": 0}

    # Rank the absolute differences with average ranks for ties.
    order = sorted(range(n), key=lambda i: abs(nz[i]))
    ranks = [0.0] * n
    i = 0
    tie_term = 0.0
    while i < n:
        j = i
        while j + 1 < n and abs(nz[order[j + 1]]) == abs(nz[order[i]]):
            j += 1
        avg = (i + 1 + j + 1) / 2.0        # average of 1-based ranks in the tie block
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        t = j - i + 1
        if t > 1:
            tie_term += t ** 3 - t
        i = j + 1

    w_plus = sum(ranks[i] for i in range(n) if nz[i] > 0)
    mean_w = n * (n + 1) / 4.0
    var_w = n * (n + 1) * (2 * n + 1) / 24.0 - tie_term / 48.0
    if var_w <= 0:
        z = 0.0
    else:
        # Continuity correction toward the mean.
        diff = w_plus - mean_w
        cc = 0.5 if diff > 0 else (-0.5 if diff < 0 else 0.0)
        z = (diff - cc) / math.sqrt(var_w)
    p = _erfc(abs(z) / math.sqrt(2.0))
    return {"statistic": w_plus, "z": z, "p_value": p, "n": n}


def sign_test(x, mu0=0.0, y=None):
    """Sign test that the (paired) sample median equals ``mu0``.

    Counts observations above ``mu0`` (or ``x_i > y_i`` for paired data), dropping
    exact ties, and tests against Binomial(n, 1/2). Returns a dict with ``n_plus``,
    ``n`` (non-tied count), the exact two-sided binomial ``p_value`` and the
    normal-approximation ``z``.
    """
    if y is not None:
        if len(x) != len(y):
            raise ValueError("x and y must have equal length")
        diffs = [x[i] - y[i] for i in range(len(x))]
    else:
        diffs = [xi - mu0 for xi in x]

    n_plus = sum(1 for d in diffs if d > 0)
    n_minus = sum(1 for d in diffs if d < 0)
    n = n_plus + n_minus
    if n == 0:
        return {"n_plus": 0, "n": 0, "p_value": 1.0, "z": 0.0}

    # Exact two-sided binomial p-value: 2 * P(K <= min(n_plus, n_minus)).
    k = min(n_plus, n_minus)
    p_two = 2.0 * binomial_cdf(k, n, 0.5)
    p_two = min(p_two, 1.0)
    # Normal approximation with continuity correction.
    mean = n / 2.0
    sd = math.sqrt(n / 4.0)
    diff = n_plus - mean
    cc = 0.5 if diff > 0 else (-0.5 if diff < 0 else 0.0)
    z = (diff - cc) / sd if sd > 0 else 0.0
    return {"n_plus": n_plus, "n": n, "p_value": p_two, "z": z}
