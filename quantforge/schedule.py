"""Payment-schedule generation with roll and business-day conventions.

Dates are ``(year, month, day)`` tuples. Builds a sequence of period end dates
from a start date, tenor in months, and total maturity in years, then applies a
business-day adjustment against a weekend calendar (Saturday/Sunday). No third-
party dependencies; uses the standard-library ``date`` for arithmetic only.
"""

from datetime import date, timedelta


def _add_months(d, n, end_of_month=False):
    """Add ``n`` months to ``(y, m, d)``, clamping the day to the month length.

    With ``end_of_month`` the result snaps to the last day of the target month
    whenever the start date is itself a month end.
    """
    y, m, day = d
    total = (m - 1) + n
    ny = y + total // 12
    nm = total % 12 + 1
    # Days in the target month.
    if nm == 12:
        last = 31
    else:
        last = (date(ny, nm + 1, 1) - timedelta(days=1)).day
    start_last = (date(y, m + 1, 1) - timedelta(days=1)).day if m < 12 else 31
    if end_of_month and day == start_last:
        nd = last
    else:
        nd = min(day, last)
    return (ny, nm, nd)


def _is_weekend(d):
    return date(d[0], d[1], d[2]).weekday() >= 5   # 5=Sat, 6=Sun


def adjust_business_day(d, convention="following"):
    """Adjust a date off weekends per a business-day convention.

    ``following`` rolls forward to the next weekday; ``preceding`` rolls back;
    ``modified_following`` rolls forward unless that crosses into the next month,
    in which case it rolls back. ``unadjusted`` returns the date unchanged.
    """
    conv = convention.lower().replace("-", "_")
    if conv == "unadjusted":
        return d
    dt = date(d[0], d[1], d[2])
    if conv in ("following", "modified_following"):
        cur = dt
        while cur.weekday() >= 5:
            cur += timedelta(days=1)
        if conv == "modified_following" and cur.month != dt.month:
            cur = dt
            while cur.weekday() >= 5:
                cur -= timedelta(days=1)
        return (cur.year, cur.month, cur.day)
    if conv == "preceding":
        cur = dt
        while cur.weekday() >= 5:
            cur -= timedelta(days=1)
        return (cur.year, cur.month, cur.day)
    raise ValueError(f"unknown business-day convention: {convention}")


def generate_schedule(start, maturity_years, freq_months, end_of_month=False,
                      convention="unadjusted"):
    """Generate period end dates from ``start`` over ``maturity_years``.

    Steps ``freq_months`` at a time (1=monthly, 3=quarterly, 6=semiannual,
    12=annual) until ``maturity_years`` is reached, then applies
    :func:`adjust_business_day` with ``convention``. Returns the list of adjusted
    period end dates (the start date itself is not included).
    """
    if freq_months < 1:
        raise ValueError("freq_months must be >= 1")
    if maturity_years <= 0:
        raise ValueError("maturity_years must be positive")
    n = int(round(maturity_years * 12 / freq_months))
    if n < 1:
        raise ValueError("schedule has no periods")
    out = []
    for i in range(1, n + 1):
        raw = _add_months(start, i * freq_months, end_of_month)
        out.append(adjust_business_day(raw, convention))
    return out
