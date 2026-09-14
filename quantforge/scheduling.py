"""Interval scheduling: maximum-weight selection, activity selection, room partition.

Three classic problems over intervals ``(start, end)`` (half-open, so ``[s, e)`` and
``[e, f)`` do not overlap):

  * ``weighted_interval_schedule`` -- pick a non-overlapping subset of maximum total weight,
    by DP after sorting on end time with a binary search for the last compatible interval
    (``O(n log n)``).
  * ``activity_selection`` -- the unweighted case: the maximum *number* of non-overlapping
    intervals, by the earliest-finish-time greedy.
  * ``min_rooms`` -- the fewest resources needed to run all intervals at once (the peak
    number simultaneously active), by an endpoint sweep.

Pure standard library.
"""

import bisect


def weighted_interval_schedule(intervals):
    """Maximum-weight non-overlapping subset. Returns ``(total_weight, chosen)``.

    ``intervals`` is a list of ``(start, end, weight)`` with ``start < end`` and
    ``weight >= 0``. ``chosen`` is the selected intervals in start order. Half-open
    intervals: one ending exactly when another starts do not conflict.
    """
    if not intervals:
        return 0.0, []
    items = sorted(intervals, key=lambda iv: iv[1])   # by end time
    ends = [iv[1] for iv in items]
    n = len(items)
    # p[i] = index of the last interval (in sorted order) that ends <= items[i].start
    best = [0.0] * (n + 1)
    take = [False] * n
    prev = [0] * n
    for i in range(n):
        s, e, w = items[i]
        if w < 0:
            raise ValueError("weights must be non-negative")
        j = bisect.bisect_right(ends, s) - 1        # last interval ending <= s
        prev[i] = j
        incl = w + best[j + 1]
        if incl > best[i]:
            best[i + 1] = incl
            take[i] = True
        else:
            best[i + 1] = best[i]
    # reconstruct
    chosen = []
    i = n - 1
    while i >= 0:
        if take[i] and best[i + 1] != best[i]:
            chosen.append(items[i])
            i = prev[i]
        else:
            i -= 1
    chosen.reverse()
    return best[n], chosen


def activity_selection(intervals):
    """Maximum number of mutually non-overlapping intervals (earliest-finish greedy).

    Returns the chosen ``(start, end)`` list in order. Ignores weights; half-open overlap.
    """
    if not intervals:
        return []
    items = sorted(intervals, key=lambda iv: (iv[1], iv[0]))
    chosen = []
    last_end = None
    for iv in items:
        s, e = iv[0], iv[1]
        if last_end is None or s >= last_end:
            chosen.append((s, e))
            last_end = e
    return chosen


def min_rooms(intervals):
    """Minimum resources to run every interval, i.e. the peak simultaneous count.

    An interval ending exactly when another begins does not need a separate room
    (half-open). Returns an integer.
    """
    if not intervals:
        return 0
    events = []
    for iv in intervals:
        s, e = iv[0], iv[1]
        events.append((s, 1))    # start: +1
        events.append((e, -1))   # end: -1 (processed before starts at the same time)
    # sort by time; at equal time, ends (-1) come before starts (+1) for half-open
    events.sort(key=lambda x: (x[0], x[1]))
    current = 0
    peak = 0
    for _, delta in events:
        current += delta
        if current > peak:
            peak = current
    return peak
