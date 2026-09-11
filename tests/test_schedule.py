"""Payment-schedule generation and business-day adjustment (schedule module)."""

import pytest

from quantforge import generate_schedule, adjust_business_day


def test_quarterly_end_of_month():
    sch = generate_schedule((2024, 1, 31), 1.0, 3, end_of_month=True)
    assert sch == [(2024, 4, 30), (2024, 7, 31), (2024, 10, 31), (2025, 1, 31)]


def test_semiannual_five_years_has_ten_periods():
    assert len(generate_schedule((2024, 1, 15), 5.0, 6)) == 10


def test_annual_last_date_is_maturity():
    sch = generate_schedule((2024, 1, 15), 3.0, 12)
    assert sch[-1] == (2027, 1, 15)
    assert len(sch) == 3


def test_month_end_clamps_short_month():
    # Jan 31 + 1 month lands on Feb 29 in a leap year.
    sch = generate_schedule((2024, 1, 31), 1.0 / 6, 1, end_of_month=True)
    assert sch[0] == (2024, 2, 29)


def test_following_convention():
    # 2024-06-01 is a Saturday -> next weekday Monday 6/3.
    assert adjust_business_day((2024, 6, 1), "following") == (2024, 6, 3)


def test_preceding_convention():
    assert adjust_business_day((2024, 6, 1), "preceding") == (2024, 5, 31)


def test_modified_following_rolls_back_across_month():
    # 2024-08-31 is a Saturday; following would cross into September, so roll back.
    assert adjust_business_day((2024, 8, 31), "modified_following") == (2024, 8, 30)


def test_unadjusted_is_identity():
    assert adjust_business_day((2024, 6, 1), "unadjusted") == (2024, 6, 1)


def test_validation():
    with pytest.raises(ValueError):
        generate_schedule((2024, 1, 1), 1.0, 0)
    with pytest.raises(ValueError):
        generate_schedule((2024, 1, 1), -1.0, 3)
    with pytest.raises(ValueError):
        adjust_business_day((2024, 6, 1), "bogus")
