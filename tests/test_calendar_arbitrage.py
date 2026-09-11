"""Model-free calendar-arbitrage check across expiries (calendar_arbitrage_violations)."""

import pytest

from quantforge import (
    calendar_arbitrage_violations, surface_is_calendar_arbitrage_free,
)
from quantforge.ssvi import SSVIParams, ssvi_calendar_free


TS = [0.25, 0.5, 1.0]


def test_flat_vol_term_structure_is_free():
    # Constant vol -> total variance w = sigma^2 t rises with t: no calendar arb.
    fns = [lambda k: 0.2, lambda k: 0.2, lambda k: 0.2]
    assert surface_is_calendar_arbitrage_free(TS, fns)
    assert calendar_arbitrage_violations(TS, fns) == []


def test_dropping_total_variance_is_flagged():
    # High short-dated vol then low long-dated vol makes w(t2) < w(t1): arbitrage.
    fns = [lambda k: 0.4, lambda k: 0.2, lambda k: 0.2]
    assert not surface_is_calendar_arbitrage_free(TS, fns)
    assert len(calendar_arbitrage_violations(TS, fns)) > 0


def test_agrees_with_ssvi_calendar_check():
    p = SSVIParams(rho=-0.3, eta=0.5, gamma=0.4,
                   thetas={0.25: 0.05, 0.5: 0.1, 1.0: 0.2})
    fns = [(lambda k, tt=tt: p.implied_vol(k, tt)) for tt in TS]
    assert surface_is_calendar_arbitrage_free(TS, fns) == ssvi_calendar_free(p)
    assert surface_is_calendar_arbitrage_free(TS, fns)


def test_violation_tuples_are_within_grid():
    fns = [lambda k: 0.4, lambda k: 0.2, lambda k: 0.2]
    for k, t_lo, t_hi in calendar_arbitrage_violations(TS, fns):
        assert -1.5 - 1e-9 <= k <= 1.5 + 1e-9
        assert t_hi > t_lo


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        calendar_arbitrage_violations([0.25, 0.5], [lambda k: 0.2])


def test_non_increasing_expiries_raise():
    with pytest.raises(ValueError):
        calendar_arbitrage_violations([0.5, 0.25], [lambda k: 0.2, lambda k: 0.2])
