"""Coin change and subset sum: classic dynamic-programming counting and reachability.

Three related DP problems over a set of positive integer denominations and a target:

  * ``min_coins`` -- the fewest coins (with unlimited supply of each) that sum to the
    target, and one such multiset.
  * ``count_change`` -- how many *distinct* multisets of coins sum to the target (order
    ignored).
  * ``subset_sum`` -- whether a subset of a given list (each item used at most once) sums to
    the target, and one witnessing subset.

All run in ``O(target * n)`` time. Pure standard library.
"""


def min_coins(coins, target):
    """Fewest coins summing to ``target`` with unlimited supply. Returns ``(count, multiset)``.

    ``coins`` are positive integer denominations. Returns ``(-1, [])`` if the target cannot
    be made. ``target == 0`` returns ``(0, [])``. The multiset is sorted ascending.
    """
    if target < 0:
        raise ValueError("target must be non-negative")
    denoms = [c for c in coins]
    for c in denoms:
        if c <= 0:
            raise ValueError("coin denominations must be positive")
    INF = float("inf")
    best = [0] + [INF] * target
    pick = [-1] * (target + 1)          # last coin used to reach each amount
    for amount in range(1, target + 1):
        for c in denoms:
            if c <= amount and best[amount - c] + 1 < best[amount]:
                best[amount] = best[amount - c] + 1
                pick[amount] = c
    if best[target] == INF:
        return -1, []
    multiset = []
    a = target
    while a > 0:
        c = pick[a]
        multiset.append(c)
        a -= c
    multiset.sort()
    return int(best[target]), multiset


def count_change(coins, target):
    """Number of distinct multisets of ``coins`` (unlimited supply) summing to ``target``.

    Order does not matter (``1 + 2`` and ``2 + 1`` count once). ``target == 0`` returns 1
    (the empty multiset).
    """
    if target < 0:
        raise ValueError("target must be non-negative")
    denoms = [c for c in coins]
    for c in denoms:
        if c <= 0:
            raise ValueError("coin denominations must be positive")
    ways = [1] + [0] * target
    for c in denoms:                    # coins in the outer loop -> multisets, not sequences
        for amount in range(c, target + 1):
            ways[amount] += ways[amount - c]
    return ways[target]


def subset_sum(values, target):
    """Whether a subset of ``values`` (each used once) sums to ``target``. Returns ``(bool, subset)``.

    ``values`` are non-negative integers. On success the second element is one witnessing
    subset (as a list of the chosen values); on failure it is ``[]``. ``target == 0`` is
    always reachable by the empty subset.
    """
    if target < 0:
        raise ValueError("target must be non-negative")
    vals = list(values)
    for v in vals:
        if v < 0:
            raise ValueError("values must be non-negative")
    # reachable[a] holds the index of the item that first completed amount a (or None)
    reachable = [False] * (target + 1)
    reachable[0] = True
    came_from = [None] * (target + 1)   # (prev_amount, value) used to reach this amount
    for idx, v in enumerate(vals):
        if v == 0:
            continue
        # iterate downward so each item is used at most once
        for a in range(target, v - 1, -1):
            if not reachable[a] and reachable[a - v]:
                reachable[a] = True
                came_from[a] = (a - v, v)
    if not reachable[target]:
        return False, []
    subset = []
    a = target
    while a > 0:
        prev, v = came_from[a]
        subset.append(v)
        a = prev
    subset.sort()
    return True, subset
