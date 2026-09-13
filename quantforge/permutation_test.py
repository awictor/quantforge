"""General permutation tests: assumption-free significance for any statistic.

A permutation test asks: if the group labels were meaningless, how often would chance
produce a group difference as extreme as the observed one? It shuffles the pooled data,
recomputes the statistic each time, and reads the p-value off that empirical null -- no
distributional assumption, valid for any test statistic (mean difference, median
difference, correlation, anything). Two forms:

  * ``permutation_test`` -- two independent samples; reshuffle the group labels.
  * ``paired_permutation_test`` -- matched pairs; flip the sign of each difference
    (equivalent to swapping the pair's two measurements).

Deterministic via a seeded LCG. Pure standard library.
"""


def _lcg(seed):
    state = seed & 0x7FFFFFFF
    def rand():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state / 0x80000000
    return rand


def _mean_diff(a, b):
    return sum(a) / len(a) - sum(b) / len(b)


def permutation_test(a, b, statistic=None, n_permutations=9999,
                     alternative="two-sided", seed=1234567):
    """Two-sample permutation test.

    ``statistic(a, b)`` defaults to the difference in means. Pools the two samples,
    reshuffles the labels ``n_permutations`` times, and returns a dict with the
    ``observed`` statistic and the ``p_value``. ``alternative`` is ``"two-sided"``,
    ``"greater"`` or ``"less"``. Deterministic for a fixed ``seed``.
    """
    if statistic is None:
        statistic = _mean_diff
    na, nb = len(a), len(b)
    if na == 0 or nb == 0:
        raise ValueError("both samples must be non-empty")
    observed = statistic(list(a), list(b))
    pool = list(a) + list(b)
    n = na + nb
    rand = _lcg(seed)

    count = 0
    for _ in range(n_permutations):
        p = pool[:]
        for i in range(n - 1, 0, -1):
            j = int(rand() * (i + 1))
            p[i], p[j] = p[j], p[i]
        stat = statistic(p[:na], p[na:])
        if alternative == "two-sided":
            if abs(stat) >= abs(observed) - 1e-15:
                count += 1
        elif alternative == "greater":
            if stat >= observed - 1e-15:
                count += 1
        elif alternative == "less":
            if stat <= observed + 1e-15:
                count += 1
        else:
            raise ValueError("alternative must be two-sided/greater/less")
    p_value = (1.0 + count) / (1.0 + n_permutations)
    return {"observed": observed, "p_value": p_value, "n_permutations": n_permutations}


def paired_permutation_test(x, y, n_permutations=9999, alternative="two-sided",
                            seed=1234567):
    """Paired permutation test on the within-pair differences ``x_i - y_i``.

    Under the null the sign of each difference is exchangeable, so each permutation
    flips signs at random. The statistic is the mean difference. Returns the same dict
    shape as :func:`permutation_test`.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n == 0:
        raise ValueError("need at least one pair")
    d = [x[i] - y[i] for i in range(n)]
    observed = sum(d) / n
    rand = _lcg(seed)

    count = 0
    for _ in range(n_permutations):
        s = sum((di if rand() < 0.5 else -di) for di in d) / n
        if alternative == "two-sided":
            if abs(s) >= abs(observed) - 1e-15:
                count += 1
        elif alternative == "greater":
            if s >= observed - 1e-15:
                count += 1
        elif alternative == "less":
            if s <= observed + 1e-15:
                count += 1
        else:
            raise ValueError("alternative must be two-sided/greater/less")
    p_value = (1.0 + count) / (1.0 + n_permutations)
    return {"observed": observed, "p_value": p_value, "n_permutations": n_permutations}
