"""Claims-triangle conversion utilities."""

from quantforge import (
    incremental_to_cumulative, cumulative_to_incremental, paid_to_date,
)


INC = [[100, 50, 30], [110, 55], [120]]
CUM = [[100, 150, 180], [110, 165], [120]]


def test_incremental_to_cumulative():
    assert incremental_to_cumulative(INC) == [[100.0, 150.0, 180.0],
                                              [110.0, 165.0], [120.0]]


def test_round_trip_both_ways():
    cum = incremental_to_cumulative(INC)
    assert cumulative_to_incremental(cum) == [[100.0, 50.0, 30.0],
                                              [110.0, 55.0], [120.0]]
    assert incremental_to_cumulative(cumulative_to_incremental(cum)) == cum


def test_cumulative_is_monotone():
    cum = incremental_to_cumulative(INC)
    assert all(all(row[j] <= row[j + 1] for j in range(len(row) - 1)) for row in cum)


def test_paid_to_date_is_last_column():
    assert paid_to_date(CUM) == [180, 165, 120]


def test_increment_sum_equals_cumulative_last():
    cum = incremental_to_cumulative(INC)
    assert all(abs(sum(INC[i]) - cum[i][-1]) < 1e-9 for i in range(len(INC)))
