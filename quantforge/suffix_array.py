"""Suffix array and LCP array: indexed substring queries over a text.

The suffix array is the permutation of a text's starting positions that lists its suffixes
in lexicographic order. Built here by prefix doubling in ``O(n log n)``, it supports
substring search by binary search in ``O(m log n)`` and, paired with the LCP (longest
common prefix) array from Kasai's ``O(n)`` algorithm, answers repeated-substring and
distinct-substring questions directly. Pure standard library.
"""


def suffix_array(text):
    """Return the suffix array of ``text`` -- suffix start positions in sorted order.

    ``result[r]`` is the starting index of the ``r``-th smallest suffix. Built by prefix
    doubling: sort by first character, then repeatedly refine ranks using pairs of ranks a
    power-of-two apart, so ``O(log n)`` rounds of ``O(n)`` counting each.
    """
    n = len(text)
    if n == 0:
        return []
    if n == 1:
        return [0]
    sa = list(range(n))
    rank = [ord(c) for c in text]
    tmp = [0] * n
    k = 1
    while True:
        # sort by (rank[i], rank[i+k]) using the pair as the key
        def key(i):
            second = rank[i + k] if i + k < n else -1
            return (rank[i], second)

        sa.sort(key=key)
        tmp[sa[0]] = 0
        for j in range(1, n):
            tmp[sa[j]] = tmp[sa[j - 1]] + (1 if key(sa[j]) != key(sa[j - 1]) else 0)
        rank = tmp[:]
        if rank[sa[n - 1]] == n - 1:      # all ranks distinct -> fully sorted
            break
        k *= 2
    return sa


def rank_array(sa):
    """Inverse of the suffix array: ``rank[i]`` is the sorted position of suffix ``i``."""
    n = len(sa)
    rank = [0] * n
    for r, i in enumerate(sa):
        rank[i] = r
    return rank


def lcp_array(text, sa=None):
    """Kasai LCP array: ``lcp[r]`` = longest common prefix of ``sa[r]`` and ``sa[r-1]``.

    ``lcp[0]`` is 0 by convention. Runs in ``O(n)`` given the suffix array (computed if not
    supplied), exploiting that adjacent suffixes in text order lose at most one leading
    character of shared prefix between successive positions.
    """
    n = len(text)
    if n == 0:
        return []
    if sa is None:
        sa = suffix_array(text)
    rank = rank_array(sa)
    lcp = [0] * n
    h = 0
    for i in range(n):
        if rank[i] > 0:
            j = sa[rank[i] - 1]          # previous suffix in sorted order
            while i + h < n and j + h < n and text[i + h] == text[j + h]:
                h += 1
            lcp[rank[i]] = h
            if h > 0:
                h -= 1
        else:
            h = 0
    return lcp


def substring_search(text, pattern, sa=None):
    """Return sorted start indices of every occurrence of ``pattern`` in ``text``.

    Binary-searches the suffix array for the block of suffixes that start with ``pattern``
    (``O(m log n)``). An empty pattern matches at every position.
    """
    n = len(text)
    m = len(pattern)
    if m == 0:
        return list(range(n + 1)) if n else [0]
    if n == 0:
        return []
    if sa is None:
        sa = suffix_array(text)

    # lower bound: first suffix >= pattern
    lo, hi = 0, n
    while lo < hi:
        mid = (lo + hi) // 2
        if text[sa[mid]:sa[mid] + m] < pattern:
            lo = mid + 1
        else:
            hi = mid
    start = lo
    # upper bound: first suffix whose prefix exceeds pattern
    hi = n
    lo = start
    while lo < hi:
        mid = (lo + hi) // 2
        if text[sa[mid]:sa[mid] + m] <= pattern:
            lo = mid + 1
        else:
            hi = mid
    end = lo
    return sorted(sa[start:end])


def longest_repeated_substring(text):
    """Return a longest substring occurring at least twice (``""`` if none repeats).

    The answer is the text slice at the maximum LCP value: the deepest shared prefix
    between two adjacent sorted suffixes. Ties resolve to the first such position.
    """
    n = len(text)
    if n < 2:
        return ""
    sa = suffix_array(text)
    lcp = lcp_array(text, sa)
    best_len = 0
    best_pos = 0
    for r in range(1, n):
        if lcp[r] > best_len:
            best_len = lcp[r]
            best_pos = sa[r]
    return text[best_pos:best_pos + best_len]


def count_distinct_substrings(text):
    """Number of distinct non-empty substrings of ``text``.

    Equals ``sum(n - sa[r]) - sum(lcp)`` -- total suffix lengths minus the prefixes shared
    with the previous sorted suffix (which would be double-counted).
    """
    n = len(text)
    if n == 0:
        return 0
    sa = suffix_array(text)
    lcp = lcp_array(text, sa)
    total = sum(n - sa[r] for r in range(n))
    return total - sum(lcp)
