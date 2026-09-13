"""Theil's U1 and U2 forecast statistics."""

import pytest

from quantforge import theil_u1, theil_u2

ACTUAL = [10, 11, 13, 12, 14, 15, 14, 16]


def test_perfect_forecast_zero():
    assert theil_u1(ACTUAL, ACTUAL) == 0.0
    assert theil_u2(ACTUAL, ACTUAL) == 0.0


def test_naive_forecast_u2_is_one():
    # Persistence forecast (each = previous actual) is the U2 benchmark.
    fc = [ACTUAL[0]] + [ACTUAL[i - 1] for i in range(1, len(ACTUAL))]
    assert abs(theil_u2(ACTUAL, fc, last_actual=ACTUAL[0]) - 1.0) < 1e-12


def test_good_forecast_beats_naive():
    good = [a + 0.1 for a in ACTUAL]
    assert theil_u2(ACTUAL, good) < 1.0


def test_bad_forecast_worse_than_naive():
    bad = [ACTUAL[i] + (5 if i % 2 else -5) for i in range(len(ACTUAL))]
    assert theil_u2(ACTUAL, bad) > 1.0


def test_u1_bounded():
    bad = [100, -100, 50, -50, 80, -80, 60, -60]
    u1 = theil_u1(ACTUAL, bad)
    assert 0.0 <= u1 <= 1.0


def test_u2_without_last_actual_skips_first():
    # Dropping last_actual evaluates from the second point onward.
    fc = ACTUAL[:]
    assert theil_u2(ACTUAL, fc) == 0.0


def test_validation():
    with pytest.raises(ValueError):
        theil_u1([1.0], [1.0, 2.0])
    with pytest.raises(ValueError):
        theil_u2([5.0, 5.0, 5.0], [5.0, 5.0, 5.0])   # naive zero error
