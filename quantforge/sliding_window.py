"""Sliding-window minimum and maximum via a monotonic deque.

For every window of ``k`` consecutive elements, the min (or max) is found in amortized
``O(1)`` -- ``O(n)`` overall -- by keeping a deque of indices whose values are monotonic:
candidates that can never win are discarded as they are dominated. This is the standard
linear-time answer to "rolling min/max", far faster than the ``O(n k)`` per-window scan.
Pure standard library.
"""

from collections import deque


def _sliding(values, k, want_max):
    n = len(values)
    if k < 1:
        raise ValueError("window size must be at least 1")
    if k > n:
        raise ValueError("window size must not exceed the sequence length")
    dq = deque()          # indices, values monotonic (decreasing for max, increasing for min)
    out = []
    for i, v in enumerate(values):
        # drop dominated tail candidates
        while dq and ((values[dq[-1]] <= v) if want_max else (values[dq[-1]] >= v)):
            dq.pop()
        dq.append(i)
        # drop the front if it has left the window
        if dq[0] <= i - k:
            dq.popleft()
        if i >= k - 1:
            out.append(values[dq[0]])
    return out


def sliding_window_max(values, k):
    """Maximum of each length-``k`` window; returns ``n - k + 1`` values."""
    return _sliding(list(values), k, want_max=True)


def sliding_window_min(values, k):
    """Minimum of each length-``k`` window; returns ``n - k + 1`` values."""
    return _sliding(list(values), k, want_max=False)


def sliding_window_sum(values, k):
    """Sum of each length-``k`` window in ``O(n)`` (running sum), ``n - k + 1`` values."""
    vals = list(values)
    n = len(vals)
    if k < 1:
        raise ValueError("window size must be at least 1")
    if k > n:
        raise ValueError("window size must not exceed the sequence length")
    s = sum(vals[:k])
    out = [s]
    for i in range(k, n):
        s += vals[i] - vals[i - k]
        out.append(s)
    return out
