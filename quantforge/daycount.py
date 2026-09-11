"""Day-count conventions: year fractions between two dates.

Dates are ``(year, month, day)`` tuples. Supports the common market
conventions:

  * ``act/360``   -- actual days / 360 (money-market standard)
  * ``act/365``   -- actual days / 365 (fixed)
  * ``30/360``    -- US (bond basis) 30/360
  * ``30E/360``   -- European 30/360
  * ``act/act``   -- ISDA actual/actual (splits across year boundaries)

Actual day counts use the proleptic Gregorian ordinal, so no ``datetime``
import is needed beyond the standard library helper.
"""

from datetime import date


def _ordinal(d):
    return date(d[0], d[1], d[2]).toordinal()


def _actual_days(start, end):
    return _ordinal(end) - _ordinal(start)


def _is_leap(y):
    return y % 4 == 0 and (y % 100 != 0 or y % 400 == 0)


def year_fraction(start, end, convention="act/365") -> float:
    """Year fraction between ``start`` and ``end`` under ``convention``.

    Both dates are ``(year, month, day)`` tuples with ``end >= start``. Negative
    intervals raise.
    """
    conv = convention.lower().replace("actual", "act")
    if _ordinal(end) < _ordinal(start):
        raise ValueError("end must not precede start")

    if conv in ("act/360", "act/360"):
        return _actual_days(start, end) / 360.0
    if conv in ("act/365", "act/365f", "act/365fixed"):
        return _actual_days(start, end) / 365.0
    if conv in ("30/360", "30/360us", "bond"):
        d1 = min(start[2], 30)
        d2 = end[2]
        if d1 == 30 and d2 == 31:
            d2 = 30
        return ((end[0] - start[0]) * 360
                + (end[1] - start[1]) * 30
                + (d2 - d1)) / 360.0
    if conv in ("30e/360", "30/360e", "eurobond"):
        d1 = min(start[2], 30)
        d2 = min(end[2], 30)
        return ((end[0] - start[0]) * 360
                + (end[1] - start[1]) * 30
                + (d2 - d1)) / 360.0
    if conv in ("act/act", "act/actisda", "act/act isda"):
        y1, y2 = start[0], end[0]
        if y1 == y2:
            denom = 366.0 if _is_leap(y1) else 365.0
            return _actual_days(start, end) / denom
        # Split: stub in the start year, whole middle years, stub in the end year.
        frac = 0.0
        start_year_end = (y1 + 1, 1, 1)
        frac += _actual_days(start, start_year_end) / (366.0 if _is_leap(y1) else 365.0)
        for y in range(y1 + 1, y2):
            frac += 1.0
        end_year_start = (y2, 1, 1)
        frac += _actual_days(end_year_start, end) / (366.0 if _is_leap(y2) else 365.0)
        return frac
    raise ValueError(f"unknown day-count convention: {convention}")


def day_count(start, end) -> int:
    """Actual number of days between two ``(y, m, d)`` dates."""
    if _ordinal(end) < _ordinal(start):
        raise ValueError("end must not precede start")
    return _actual_days(start, end)
