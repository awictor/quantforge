"""Effective/realized/impact spread decomposition."""

import pytest

from quantforge import (effective_spread, realized_spread, price_impact,
                        quoted_spread)

PRICES = [100.5, 99.6, 100.2, 101.0, 99.4]
MIDS = [100.0, 100.0, 100.0, 100.5, 100.0]
FUT = [100.2, 99.9, 100.1, 100.8, 99.7]
SIGNS = [1, -1, 1, 1, -1]


def test_decomposition_identity():
    es = effective_spread(PRICES, MIDS, SIGNS)
    rs = realized_spread(PRICES, MIDS, FUT, SIGNS)
    pi = price_impact(PRICES, MIDS, FUT, SIGNS)
    assert abs(es - (rs + pi)) < 1e-12


def test_effective_positive_when_trades_cross_spread():
    # Buys above mid, sells below mid -> positive effective spread.
    es = effective_spread([100.5, 99.5], [100.0, 100.0], [1, -1])
    assert es > 0


def test_quoted_spread_value():
    assert abs(quoted_spread([99.9, 99.8], [100.1, 100.2]) - 0.003) < 1e-9


def test_no_impact_makes_realized_equal_effective():
    # If the midpoint does not move, realized == effective and impact == 0.
    es = effective_spread(PRICES, MIDS, SIGNS)
    rs = realized_spread(PRICES, MIDS, MIDS, SIGNS)   # future mid == mid
    pi = price_impact(PRICES, MIDS, MIDS, SIGNS)
    assert abs(rs - es) < 1e-12
    assert abs(pi) < 1e-12


def test_permanent_move_is_pure_impact():
    # Trade exactly at the mid, but the mid moves in the trade direction:
    # effective is zero, so realized is negative and impact positive.
    es = effective_spread([100.0], [100.0], [1])
    pi = price_impact([100.0], [100.0], [100.2], [1])
    rs = realized_spread([100.0], [100.0], [100.2], [1])
    assert abs(es) < 1e-12
    assert pi > 0
    assert abs(es - (rs + pi)) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        effective_spread([100.0], [100.0], [2])          # bad sign
    with pytest.raises(ValueError):
        effective_spread([100.0], [0.0], [1])            # non-positive mid
    with pytest.raises(ValueError):
        realized_spread([100.0], [100.0], [100.0, 100.0], [1])
    with pytest.raises(ValueError):
        quoted_spread([100.0], [])
