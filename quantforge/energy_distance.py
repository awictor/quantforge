"""Energy distance: a distribution-free two-sample test of equal distributions.

Given samples ``X`` (size n) and ``Y`` (size m), the energy distance is

    E(X, Y) = 2 A - B - C,

where ``A`` is the mean pairwise ``|x_i - y_j|`` across the two samples, ``B`` the
mean within-``X`` distance, and ``C`` the mean within-``Y`` distance. It is zero if
and only if ``X`` and ``Y`` have the same distribution and strictly positive
otherwise -- unlike a t-test it sees differences in shape, spread, and tails, not
only the mean. Significance comes from a permutation test: pool the samples, reshuffle
the group labels many times, and compare. Univariate samples; pure standard library.
"""

import math


def _mean_cross(a, b):
    n, m = len(a), len(b)
    total = 0.0
    for x in a:
        for y in b:
            total += abs(x - y)
    return total / (n * m)


def _mean_within(a):
    n = len(a)
    if n < 2:
        return 0.0
    total = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            total += abs(a[i] - a[j])
    return 2.0 * total / (n * n)   # divide by n^2 (includes the zero diagonal)


def energy_distance(a, b):
    """Energy distance ``2 A - B - C`` between two samples.

    Non-negative; zero only when the empirical distributions coincide. ``A`` is the
    mean cross-sample absolute distance, ``B`` and ``C`` the mean within-sample
    distances (each normalized by the squared sample size, i.e. including the zero
    diagonal, per Szekely-Rizzo).
    """
    if len(a) == 0 or len(b) == 0:
        raise ValueError("both samples must be non-empty")
    A = _mean_cross(a, b)
    B = _mean_within(a)
    C = _mean_within(b)
    e = 2.0 * A - B - C
    return e if e > 0.0 else 0.0


def _lcg(seed):
    state = seed & 0x7FFFFFFF
    def rand():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state
    return rand


def energy_test(a, b, n_permutations=999, seed=1234567):
    """Permutation test of equal distributions via the energy distance.

    Pools the two samples, reshuffles the group labels ``n_permutations`` times, and
    returns a dict with the observed ``statistic``, the ``p_value``
    ``(1 + #{perm >= observed}) / (1 + n_permutations)``, and ``n_permutations``.
    A small p-value is evidence the two samples come from different distributions.
    """
    n = len(a)
    if n == 0 or len(b) == 0:
        raise ValueError("both samples must be non-empty")
    observed = energy_distance(a, b)
    pool = list(a) + list(b)
    total = len(pool)
    rand = _lcg(seed)

    count = 0
    for _ in range(n_permutations):
        # Fisher-Yates shuffle of a copy, then split at n.
        p = pool[:]
        for i in range(total - 1, 0, -1):
            j = rand() % (i + 1)
            p[i], p[j] = p[j], p[i]
        stat = energy_distance(p[:n], p[n:])
        if stat >= observed - 1e-15:
            count += 1
    p_value = (1.0 + count) / (1.0 + n_permutations)
    return {"statistic": observed, "p_value": p_value, "n_permutations": n_permutations}
