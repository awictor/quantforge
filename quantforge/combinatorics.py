"""Exact integer combinatorics: binomials, Stirling, Catalan, partitions, Bell.

Exact counting functions on Python big integers -- no floating error, no overflow. The
binomial and multinomial coefficients, Stirling numbers of the second kind (set
partitions into a fixed number of blocks), Catalan numbers, the integer-partition count,
Bell numbers (all set partitions), and derangements. These are the building blocks of
enumeration, probability on finite structures, and series coefficients. Pure standard
library.
"""

from math import comb, factorial


def binomial(n, k):
    """Binomial coefficient ``C(n, k)`` -- the number of ``k``-subsets of ``n`` items.

    Zero when ``k < 0`` or ``k > n``. Exact for all non-negative ``n`` (big integers).
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    if k < 0 or k > n:
        return 0
    return comb(n, k)


def multinomial(counts):
    """Multinomial coefficient ``(sum counts)! / prod(counts!)``.

    The number of distinct arrangements of a multiset with the given group ``counts``
    (e.g. ``multinomial([1, 4, 4]) == 630`` for the letters of "mississippi" is built
    from these). All counts must be non-negative.
    """
    total = 0
    for c in counts:
        if c < 0:
            raise ValueError("counts must be non-negative")
        total += c
    result = factorial(total)
    for c in counts:
        result //= factorial(c)
    return result


def stirling_second(n, k):
    """Stirling number of the second kind ``S(n, k)``: partitions of ``n`` items into ``k`` non-empty blocks.

    Uses the recurrence ``S(n, k) = k*S(n-1, k) + S(n-1, k-1)``. ``S(0, 0) = 1``.
    """
    if n < 0 or k < 0:
        raise ValueError("n and k must be >= 0")
    if k == 0:
        return 1 if n == 0 else 0
    if k > n:
        return 0
    # Row-by-row DP.
    prev = [0] * (k + 1)
    prev[0] = 1
    for i in range(1, n + 1):
        cur = [0] * (k + 1)
        for j in range(1, min(i, k) + 1):
            cur[j] = j * prev[j] + prev[j - 1]
        prev = cur
    return prev[k]


def bell(n):
    """Bell number ``B(n)``: the total number of partitions of an ``n``-element set.

    Computed by summing Stirling numbers of the second kind, ``B(n) = sum_k S(n, k)``,
    via the Bell triangle. ``B(0) = 1``.
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    if n == 0:
        return 1
    row = [1]
    for _ in range(n):
        nxt = [row[-1]]
        for v in row:
            nxt.append(nxt[-1] + v)
        row = nxt
    return row[0]


def catalan(n):
    """Catalan number ``C_n = binomial(2n, n) / (n + 1)``.

    Counts balanced-parenthesis strings, binary trees, monotone lattice paths, and many
    other structures. Exact.
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    return comb(2 * n, n) // (n + 1)


def partition_count(n):
    """Number of integer partitions of ``n`` (unordered sums of positive integers).

    ``p(0) = 1``; ``p(4) = 5`` (4, 3+1, 2+2, 2+1+1, 1+1+1+1). Uses the standard
    coin-change DP over parts ``1..n``, so it is exact and ``O(n^2)``.
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    dp = [0] * (n + 1)
    dp[0] = 1
    for part in range(1, n + 1):
        for total in range(part, n + 1):
            dp[total] += dp[total - part]
    return dp[n]


def derangements(n):
    """Number of derangements ``D(n)`` -- permutations of ``n`` items with no fixed point.

    Uses the recurrence ``D(n) = (n-1) * (D(n-1) + D(n-2))`` with ``D(0) = 1``,
    ``D(1) = 0``. Exact.
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    if n == 0:
        return 1
    if n == 1:
        return 0
    a, b = 1, 0                      # D(0), D(1)
    for i in range(2, n + 1):
        a, b = b, (i - 1) * (b + a)
    return b
