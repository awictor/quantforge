"""Tests for calendar arithmetic, cross-checked against Python's datetime/calendar."""

import calendar as pycal
import datetime
import random

import pytest

from quantforge.calendar_math import (
    is_leap_year,
    days_in_month,
    julian_day_number,
    jdn_to_date,
    days_between,
    day_of_week,
    day_of_week_name,
    day_of_year,
    add_days,
    easter_date,
)


def test_fuzz_vs_datetime():
    rng = random.Random(221)
    for _ in range(8000):
        y = rng.randint(1, 9999)
        m = rng.randint(1, 12)
        d = rng.randint(1, pycal.monthrange(y, m)[1])
        dt = datetime.date(y, m, d)
        assert day_of_week(y, m, d) == dt.weekday()
        assert day_of_year(y, m, d) == dt.timetuple().tm_yday
        assert jdn_to_date(julian_day_number(y, m, d)) == (y, m, d)
        assert is_leap_year(y) == pycal.isleap(y)
        assert days_in_month(y, m) == pycal.monthrange(y, m)[1]
        lo = (datetime.date(1, 1, 1) - dt).days
        hi = (datetime.date(9999, 12, 31) - dt).days
        n = rng.randint(max(-3000, lo), min(3000, hi))
        dt2 = dt + datetime.timedelta(days=n)
        assert add_days(y, m, d, n) == (dt2.year, dt2.month, dt2.day)
        y2 = rng.randint(1, 9999)
        m2 = rng.randint(1, 12)
        d2 = rng.randint(1, pycal.monthrange(y2, m2)[1])
        assert days_between((y, m, d), (y2, m2, d2)) == (datetime.date(y2, m2, d2) - dt).days


def test_known_easter_dates():
    known = {
        2000: (2000, 4, 23),
        2020: (2020, 4, 12),
        2024: (2024, 3, 31),
        2025: (2025, 4, 20),
        1961: (1961, 4, 2),
        1943: (1943, 4, 25),
    }
    for yr, e in known.items():
        assert easter_date(yr) == e


def test_jdn_reference():
    assert julian_day_number(2000, 1, 1) == 2451545  # standard reference epoch


def test_weekday_names():
    assert day_of_week_name(2024, 1, 1) == "Monday"
    assert day_of_week_name(2000, 1, 1) == "Saturday"
    assert day_of_week_name(2026, 9, 9) == "Wednesday"


def test_leap_year_rules():
    assert is_leap_year(2000) is True   # divisible by 400
    assert is_leap_year(1900) is False  # divisible by 100, not 400
    assert is_leap_year(2024) is True
    assert is_leap_year(2023) is False


def test_days_in_month():
    assert days_in_month(2024, 2) == 29
    assert days_in_month(2023, 2) == 28
    assert days_in_month(2024, 4) == 30
    assert days_in_month(2024, 12) == 31


def test_add_days_proleptic_beyond_datetime_range():
    # module is not bounded by datetime's year <= 9999
    assert add_days(9999, 12, 31, 1) == (10000, 1, 1)
    assert add_days(2024, 1, 1, -1) == (2023, 12, 31)


def test_days_between_signs():
    assert days_between((2024, 1, 1), (2024, 1, 11)) == 10
    assert days_between((2024, 1, 11), (2024, 1, 1)) == -10
    assert days_between((2024, 1, 1), (2024, 1, 1)) == 0


def test_day_of_year_bounds():
    assert day_of_year(2024, 1, 1) == 1
    assert day_of_year(2024, 12, 31) == 366  # leap
    assert day_of_year(2023, 12, 31) == 365


def test_round_trip_over_epochs():
    for jdn in (0, 1721426, 2451545, 2460000, 5000000):
        assert julian_day_number(*jdn_to_date(jdn)) == jdn


def test_invalid_month_raises():
    with pytest.raises(ValueError):
        julian_day_number(2020, 13, 1)
    with pytest.raises(ValueError):
        days_in_month(2020, 0)


def test_invalid_day_raises():
    with pytest.raises(ValueError):
        julian_day_number(2021, 2, 29)  # 2021 not a leap year
    with pytest.raises(ValueError):
        julian_day_number(2024, 4, 31)  # April has 30 days


def test_feb29_valid_in_leap_year():
    assert jdn_to_date(julian_day_number(2024, 2, 29)) == (2024, 2, 29)
