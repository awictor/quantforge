"""Calendar arithmetic from first principles: Julian day numbers, weekdays, Easter.

Dates in the proleptic Gregorian calendar reduce to integer arithmetic once mapped to a
*Julian day number* (JDN) -- a continuous count of days that makes differences, weekdays,
and round-trips trivial. This module implements the standard Fliegel-Van Flandern
conversion, plus leap-year testing, day-of-week and day-of-year, and the Gregorian Easter
(the "Computus"). Dates are ``(year, month, day)`` tuples. No ``datetime`` dependency. Pure
standard library.
"""

_WEEKDAYS = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
_MONTH_DAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def is_leap_year(year):
    """True if ``year`` is a Gregorian leap year (divisible by 4, not 100 unless 400)."""
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def days_in_month(year, month):
    """Number of days in ``month`` of ``year`` (1-12)."""
    if not 1 <= month <= 12:
        raise ValueError("month must be in 1..12")
    if month == 2 and is_leap_year(year):
        return 29
    return _MONTH_DAYS[month - 1]


def _validate(year, month, day):
    if not 1 <= month <= 12:
        raise ValueError("month must be in 1..12")
    if not 1 <= day <= days_in_month(year, month):
        raise ValueError("day out of range for the given month")


def julian_day_number(year, month, day):
    """Julian day number of a proleptic-Gregorian date (Fliegel-Van Flandern).

    A monotone integer count of days, so ``jdn(b) - jdn(a)`` is the day difference.
    """
    _validate(year, month, day)
    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    return (
        day
        + (153 * m + 2) // 5
        + 365 * y
        + y // 4
        - y // 100
        + y // 400
        - 32045
    )


def jdn_to_date(jdn):
    """Inverse of :func:`julian_day_number`: the ``(year, month, day)`` for a JDN."""
    a = jdn + 32044
    b = (4 * a + 3) // 146097
    c = a - (146097 * b) // 4
    d = (4 * c + 3) // 1461
    e = c - (1461 * d) // 4
    m = (5 * e + 2) // 153
    day = e - (153 * m + 2) // 5 + 1
    month = m + 3 - 12 * (m // 10)
    year = 100 * b + d - 4800 + m // 10
    return (year, month, day)


def days_between(date1, date2):
    """Signed day count ``date2 - date1`` (positive when ``date2`` is later)."""
    return julian_day_number(*date2) - julian_day_number(*date1)


def day_of_week(year, month, day):
    """Weekday index, 0 = Monday .. 6 = Sunday."""
    # JDN mod 7 gives the Monday-based index directly (JDN 2451545 = 2000-01-01, a Saturday)
    return julian_day_number(year, month, day) % 7


def day_of_week_name(year, month, day):
    """Weekday name, e.g. ``"Wednesday"``."""
    return _WEEKDAYS[day_of_week(year, month, day)]


def day_of_year(year, month, day):
    """Ordinal day within the year (Jan 1 = 1)."""
    return julian_day_number(year, month, day) - julian_day_number(year, 1, 1) + 1


def add_days(year, month, day, n):
    """Return the date ``n`` days after (or before, if ``n < 0``) the given date."""
    return jdn_to_date(julian_day_number(year, month, day) + n)


def easter_date(year):
    """Gregorian Easter Sunday of ``year`` as ``(year, month, day)`` (Anonymous Computus)."""
    a = year % 19
    b = year // 100
    c = year % 100
    d = b // 4
    e = b % 4
    f = (b + 8) // 25
    g = (b - f + 1) // 3
    h = (19 * a + b - d - g + 15) % 30
    i = c // 4
    k = c % 4
    l = (32 + 2 * e + 2 * i - h - k) % 7
    m = (a + 11 * h + 22 * l) // 451
    month = (h + l - 7 * m + 114) // 31
    day = ((h + l - 7 * m + 114) % 31) + 1
    return (year, month, day)
