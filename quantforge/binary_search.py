"""Binary and ternary search on predicates and integer ranges.

Search primitives that turn a monotone or unimodal property into an ``O(log n)`` query:

  * ``first_true`` -- the smallest integer in ``[lo, hi]`` where a monotone predicate
    (``False ... False True ... True``) first holds; the workhorse behind "binary search the
    answer".
  * ``last_true`` -- the largest integer where a ``True ... True False ... False`` predicate
    still holds.
  * ``ternary_search_int_max`` / ``ternary_search_int_min`` -- the integer that maximizes or
    minimizes a strictly unimodal function over ``[lo, hi]``.

Pure standard library.
"""


def first_true(lo, hi, predicate):
    """Smallest integer ``x`` in ``[lo, hi]`` with ``predicate(x)`` true, or ``hi + 1`` if none.

    ``predicate`` must be monotone: once true it stays true. ``O(log(hi - lo))`` calls.
    """
    if lo > hi:
        raise ValueError("lo must not exceed hi")
    result = hi + 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if predicate(mid):
            result = mid
            hi = mid - 1
        else:
            lo = mid + 1
    return result


def last_true(lo, hi, predicate):
    """Largest integer ``x`` in ``[lo, hi]`` with ``predicate(x)`` true, or ``lo - 1`` if none.

    ``predicate`` must be monotone-decreasing: once false it stays false.
    """
    if lo > hi:
        raise ValueError("lo must not exceed hi")
    result = lo - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if predicate(mid):
            result = mid
            lo = mid + 1
        else:
            hi = mid - 1
    return result


def ternary_search_int_max(lo, hi, f):
    """Integer in ``[lo, hi]`` maximizing a strictly unimodal ``f`` (increasing then decreasing).

    Narrows the range by thirds until three or fewer candidates remain, then takes the best.
    ``O(log(hi - lo))`` evaluations.
    """
    if lo > hi:
        raise ValueError("lo must not exceed hi")
    while hi - lo > 2:
        m1 = lo + (hi - lo) // 3
        m2 = hi - (hi - lo) // 3
        if f(m1) < f(m2):
            lo = m1 + 1
        else:
            hi = m2 - 1
    best = lo
    for x in range(lo + 1, hi + 1):
        if f(x) > f(best):
            best = x
    return best


def ternary_search_int_min(lo, hi, f):
    """Integer in ``[lo, hi]`` minimizing a strictly unimodal ``f`` (decreasing then increasing)."""
    if lo > hi:
        raise ValueError("lo must not exceed hi")
    while hi - lo > 2:
        m1 = lo + (hi - lo) // 3
        m2 = hi - (hi - lo) // 3
        if f(m1) > f(m2):
            lo = m1 + 1
        else:
            hi = m2 - 1
    best = lo
    for x in range(lo + 1, hi + 1):
        if f(x) < f(best):
            best = x
    return best
