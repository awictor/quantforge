"""Operations on sets of 1-D intervals: merge, union, intersection, overlap.

Given a collection of ``(start, end)`` ranges, these compute the merged (non-overlapping)
cover, the union and intersection of two collections, the total length covered, and the
maximum number of intervals overlapping at any point. The workhorse for scheduling,
coverage analysis, and range bookkeeping. Intervals are closed ``[start, end]`` with
``start <= end``. Pure standard library.
"""


def _normalize(intervals):
    out = []
    for iv in intervals:
        lo, hi = iv
        if lo > hi:
            raise ValueError("interval start must be <= end")
        out.append((lo, hi))
    return out


def merge_intervals(intervals):
    """Merge overlapping/adjacent intervals into a minimal sorted non-overlapping cover.

    ``[(1,3),(2,6),(8,10)] -> [(1,6),(8,10)]``. Touching intervals (``(1,2),(2,3)``) merge
    into ``(1,3)``. Returns a new sorted list; empty input yields ``[]``.
    """
    ivs = sorted(_normalize(intervals))
    if not ivs:
        return []
    merged = [ivs[0]]
    for lo, hi in ivs[1:]:
        last_lo, last_hi = merged[-1]
        if lo <= last_hi:
            merged[-1] = (last_lo, max(last_hi, hi))
        else:
            merged.append((lo, hi))
    return merged


def total_covered_length(intervals):
    """Total length covered by the union of the intervals (overlaps counted once)."""
    return sum(hi - lo for lo, hi in merge_intervals(intervals))


def intervals_intersection(a, b):
    """Intersection of two interval collections: the ranges covered by *both*.

    Merges each side first, then sweeps for overlaps. Returns a sorted
    non-overlapping list.
    """
    ma = merge_intervals(a)
    mb = merge_intervals(b)
    out = []
    i = j = 0
    while i < len(ma) and j < len(mb):
        lo = max(ma[i][0], mb[j][0])
        hi = min(ma[i][1], mb[j][1])
        if lo <= hi:
            out.append((lo, hi))
        if ma[i][1] < mb[j][1]:
            i += 1
        else:
            j += 1
    return out


def intervals_union(a, b):
    """Union of two interval collections as a merged non-overlapping list."""
    return merge_intervals(list(a) + list(b))


def max_overlap(intervals):
    """Maximum number of intervals overlapping at any single point.

    Sweeps the endpoints (starts before ends at the same coordinate, so touching
    intervals count as overlapping). Returns the peak count; ``0`` for no intervals.
    """
    events = []
    for lo, hi in _normalize(intervals):
        events.append((lo, 1))
        events.append((hi, -1))
    # Sort by coordinate; at a tie, process starts (+1) before ends (-1).
    events.sort(key=lambda e: (e[0], -e[1]))
    best = 0
    current = 0
    for _, delta in events:
        current += delta
        best = max(best, current)
    return best
