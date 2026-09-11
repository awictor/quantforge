"""Full surface static-arbitrage report: butterfly + calendar (surface_arbitrage_report)."""

import pytest

from quantforge import (
    surface_arbitrage_report, surface_is_arbitrage_free,
)
from quantforge.ssvi import SSVIParams, ssvi_is_arbitrage_free
from quantforge.svi import SVIParams


S0, R = 100.0, 0.02
TS = [0.25, 0.5, 1.0]


def _ssvi():
    return SSVIParams(rho=-0.3, eta=0.5, gamma=0.4,
                      thetas={0.25: 0.05, 0.5: 0.1, 1.0: 0.2})


def _ssvi_fns(p):
    return [(lambda k, tt=tt: p.implied_vol(k, tt)) for tt in TS]


def test_arbitrage_free_ssvi_surface():
    p = _ssvi()
    assert ssvi_is_arbitrage_free(p)                      # SSVI's own check
    assert surface_is_arbitrage_free(S0, R, TS, _ssvi_fns(p))
    rep = surface_arbitrage_report(S0, R, TS, _ssvi_fns(p))
    assert rep["butterfly"] == {}
    assert rep["calendar"] == []


def test_butterfly_violating_slice_reported():
    p = _ssvi()
    bad = SVIParams(a=0.005, b=0.3, rho=-0.95, m=0.0, s=0.01)
    fns = [(lambda k: bad.implied_vol(k, 0.25)),
           (lambda k: p.implied_vol(k, 0.5)),
           (lambda k: p.implied_vol(k, 1.0))]
    rep = surface_arbitrage_report(S0, R, TS, fns)
    assert 0.25 in rep["butterfly"]
    assert not surface_is_arbitrage_free(S0, R, TS, fns)


def test_calendar_violating_surface_reported():
    fns = [lambda k: 0.4, lambda k: 0.2, lambda k: 0.2]
    rep = surface_arbitrage_report(S0, R, TS, fns)
    assert len(rep["calendar"]) > 0
    assert not surface_is_arbitrage_free(S0, R, TS, fns)


def test_flat_surface_is_free():
    fns = [lambda k: 0.2, lambda k: 0.2, lambda k: 0.2]
    assert surface_is_arbitrage_free(S0, R, TS, fns)


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        surface_arbitrage_report(S0, R, TS, [lambda k: 0.2])
