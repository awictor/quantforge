"""Day-count conventions / year fractions (daycount module)."""

import pytest

from quantforge import year_fraction, day_count


def test_act360_half_year():
    assert year_fraction((2024, 1, 1), (2024, 6, 29), "act/360") == pytest.approx(
        180 / 360, abs=1e-12)
    assert day_count((2024, 1, 1), (2024, 6, 29)) == 180


def test_act365_leap_year_full():
    # 2024 is a leap year: 366 actual days over 365.
    assert year_fraction((2024, 1, 1), (2025, 1, 1), "act/365") == pytest.approx(
        366 / 365, abs=1e-9)


def test_act_act_same_leap_year_is_one():
    assert year_fraction((2024, 1, 1), (2025, 1, 1), "act/act") == pytest.approx(
        1.0, abs=1e-12)


def test_30_360_clean_year():
    assert year_fraction((2024, 1, 15), (2025, 1, 15), "30/360") == pytest.approx(
        1.0, abs=1e-12)


def test_30_360_month_end_rule():
    # Jan 31 -> Feb 28 counts as 28 days under 30/360 US.
    assert year_fraction((2024, 1, 31), (2024, 2, 28), "30/360") * 360 == pytest.approx(
        28.0, abs=1e-9)


def test_30e_360_caps_day_at_30():
    # Jan 31 -> Mar 31 both cap to 30 -> exactly 2 months.
    assert year_fraction((2024, 1, 31), (2024, 3, 31), "30e/360") == pytest.approx(
        60 / 360, abs=1e-12)


def test_act_act_spans_year_boundary():
    # Non-leap 2023 stub + leap 2024 stub, ~1 year but not exactly.
    v = year_fraction((2023, 7, 1), (2024, 7, 1), "act/act")
    assert 1.0 < v < 1.01


def test_validation():
    with pytest.raises(ValueError):
        year_fraction((2024, 6, 1), (2024, 1, 1))     # end before start
    with pytest.raises(ValueError):
        year_fraction((2024, 1, 1), (2024, 2, 1), "bogus")
    with pytest.raises(ValueError):
        day_count((2024, 6, 1), (2024, 1, 1))
