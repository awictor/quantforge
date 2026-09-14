"""Weighted sampling from a stream: Efraimidis-Spirakis reservoir and with-replacement.

Selecting ``k`` items with probability proportional to weight, in one pass over a stream
of unknown length: the Efraimidis-Spirakis A-Res algorithm assigns each item a key
``u^(1/w)`` and keeps the ``k`` largest -- weighted sampling *without* replacement in
``O(n log k)``. Also provides weighted sampling *with* replacement via the alias-free
cumulative method. Deterministic seeded stream for reproducibility. Pure standard library.
"""

import heapq


def _lcg(seed):
    state = seed & 0x7FFFFFFF or 1

    def _next():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return (state + 0.5) / 0x80000000

    return _next


def weighted_reservoir_sample(items, weights, k, seed=1234567):
    """Sample ``k`` distinct items with weight-proportional probability (Efraimidis-Spirakis).

    Assigns each item the key ``u^(1/w)`` for a uniform ``u`` and keeps the ``k`` largest
    keys -- one pass, ``O(n log k)`` memory ``O(k)``. Returns a list of the chosen items
    (order not significant). ``k`` is clamped to the number of positive-weight items.
    """
    if len(items) != len(weights):
        raise ValueError("items and weights must have equal length")
    if k < 0:
        raise ValueError("k must be non-negative")
    if k == 0:
        return []
    rand = _lcg(seed)
    heap = []                                    # min-heap of (key, index)
    for i, (it, w) in enumerate(zip(items, weights)):
        if w <= 0:
            continue
        u = rand()
        key = u ** (1.0 / w)
        if len(heap) < k:
            heapq.heappush(heap, (key, i, it))
        elif key > heap[0][0]:
            heapq.heapreplace(heap, (key, i, it))
    return [it for _, _, it in heap]


def weighted_sample_with_replacement(items, weights, k, seed=1234567):
    """Sample ``k`` items *with* replacement, weight-proportional (cumulative search).

    Each of the ``k`` draws is independent, so an item can appear multiple times. Returns
    a list of length ``k``. Weights must be non-negative and not all zero.
    """
    n = len(items)
    if n != len(weights):
        raise ValueError("items and weights must have equal length")
    if k < 0:
        raise ValueError("k must be non-negative")
    cum = []
    total = 0.0
    for w in weights:
        if w < 0:
            raise ValueError("weights must be non-negative")
        total += w
        cum.append(total)
    if total <= 0.0:
        raise ValueError("weights must sum to a positive value")
    rand = _lcg(seed)
    import bisect
    out = []
    for _ in range(k):
        r = rand() * total
        idx = bisect.bisect_right(cum, r)
        if idx >= n:
            idx = n - 1
        out.append(items[idx])
    return out
