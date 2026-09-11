"""Overnight compounding with lookback / lockout (rates module)."""

import pytest

from quantforge import (
    compounded_overnight_rate, compounded_rate_with_lookback,
    compounded_rate_with_lockout,
)


N = 10
ACC = [1 / 360] * N
FIX = [0.03 + 0.001 * i for i in range(N)]   # rising fixings


def test_zero_lookback_reduces_to_base():
    assert compounded_rate_with_lookback(FIX, ACC, 0) == pytest.approx(
        compounded_overnight_rate(FIX, ACC), abs=1e-12)


def test_zero_lockout_reduces_to_base():
    assert compounded_rate_with_lockout(FIX, ACC, 0) == pytest.approx(
        compounded_overnight_rate(FIX, ACC), abs=1e-12)


def test_lockout_lower_on_rising_fixings():
    # Freezing the last days at an earlier (lower) rate lowers the coupon.
    assert compounded_rate_with_lockout(FIX, ACC, 3) < compounded_overnight_rate(FIX, ACC)


def test_lookback_lower_on_rising_fixings():
    assert compounded_rate_with_lookback(FIX, ACC, 3) < compounded_overnight_rate(FIX, ACC)


def test_lockout_freezes_last_rate():
    # With full lockout of n-1 the whole period uses the first fixing.
    r = compounded_rate_with_lockout(FIX, ACC, N - 1)
    flat = compounded_overnight_rate([FIX[0]] * N, ACC)
    assert r == pytest.approx(flat, abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        compounded_rate_with_lookback(FIX, ACC, -1)
    with pytest.raises(ValueError):
        compounded_rate_with_lockout(FIX[:2], ACC[:2], 5)
    with pytest.raises(ValueError):
        compounded_rate_with_lookback([0.05], [1 / 360, 1 / 360], 0)
