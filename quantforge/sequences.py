"""Classic sequence algorithms: LIS, maximum subarray, longest run.

The longest increasing subsequence (patience-sorting, ``O(n log n)``), the maximum-sum
contiguous subarray (Kadane's algorithm), and the longest run of a repeated value. Each
returns both the answer and the witnessing indices/values where useful. Pure standard
library.
"""

import bisect


def longest_increasing_subsequence(x, strict=True):
    """Longest (strictly by default) increasing subsequence of ``x``.

    Returns an actual longest subsequence as a list (not necessarily contiguous). With
    ``strict=False`` allows equal consecutive values (non-decreasing). ``O(n log n)`` via
    patience sorting with predecessor tracking. Empty input yields ``[]``.
    """
    n = len(x)
    if n == 0:
        return []
    tails = []            # tails[i] = index in x of the smallest tail of an LIS of length i+1
    prev = [-1] * n
    tail_vals = []        # parallel list of the tail values for bisect
    for i, v in enumerate(x):
        if strict:
            pos = bisect.bisect_left(tail_vals, v)
        else:
            pos = bisect.bisect_right(tail_vals, v)
        if pos == len(tails):
            tails.append(i)
            tail_vals.append(v)
        else:
            tails[pos] = i
            tail_vals[pos] = v
        prev[i] = tails[pos - 1] if pos > 0 else -1
    # Reconstruct from the last tail.
    seq = []
    k = tails[-1]
    while k != -1:
        seq.append(x[k])
        k = prev[k]
    seq.reverse()
    return seq


def maximum_subarray(x):
    """Maximum-sum contiguous subarray (Kadane): ``(sum, start, end)`` inclusive.

    Returns the largest achievable sum of a non-empty contiguous slice and its bounds.
    Handles all-negative inputs (returns the single largest element). Raises on empty
    input.
    """
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    best_sum = x[0]
    best_lo = best_hi = 0
    cur = x[0]
    cur_lo = 0
    for i in range(1, n):
        if cur < 0:
            cur = x[i]
            cur_lo = i
        else:
            cur += x[i]
        if cur > best_sum:
            best_sum = cur
            best_lo = cur_lo
            best_hi = i
    return best_sum, best_lo, best_hi


def longest_run(x):
    """Longest run of a single repeated value: ``(value, length, start_index)``.

    Scans for the longest maximal streak of equal adjacent elements. Ties break to the
    earliest run. Raises on empty input.
    """
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    best_val = x[0]
    best_len = 1
    best_start = 0
    cur_len = 1
    cur_start = 0
    for i in range(1, n):
        if x[i] == x[i - 1]:
            cur_len += 1
        else:
            cur_len = 1
            cur_start = i
        if cur_len > best_len:
            best_len = cur_len
            best_val = x[i]
            best_start = cur_start
    return best_val, best_len, best_start
