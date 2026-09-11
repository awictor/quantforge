"""Box-spread implied financing rate (bsm module)."""

import pytest

from quantforge import box_spread_implied_rate
from quantforge.strategy import box_spread


K1, K2, T, SIG = 90.0, 110.0, 1.0, 0.2


def _box_value(r):
    return box_spread(100.0, K1, K2, T, r, SIG).net.market_value


def test_recovers_financing_rate():
    val = _box_value(0.05)
    assert box_spread_implied_rate(val, K1, K2, T) == pytest.approx(0.05, abs=1e-9)


def test_recovers_other_rate():
    val = _box_value(0.02)
    assert box_spread_implied_rate(val, K1, K2, T) == pytest.approx(0.02, abs=1e-9)


def test_negative_rate_above_notional():
    # A box priced above its (K2-K1) notional implies a negative rate.
    assert box_spread_implied_rate(20.5, K1, K2, T) < 0.0


def test_at_notional_zero_rate():
    assert box_spread_implied_rate(20.0, K1, K2, T) == pytest.approx(0.0, abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        box_spread_implied_rate(19.0, 110.0, 90.0, T)   # K order
    with pytest.raises(ValueError):
        box_spread_implied_rate(19.0, K1, K2, 0.0)      # t = 0
    with pytest.raises(ValueError):
        box_spread_implied_rate(-1.0, K1, K2, T)        # negative price
