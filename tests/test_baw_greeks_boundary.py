"""BAW American critical spot + greeks (baw_critical_spot, baw_american_greeks)."""

import pytest

from quantforge import (
    baw_american, baw_critical_spot, baw_american_greeks,
)
from quantforge.bsm import call_price, put_price


def test_no_dividend_call_has_no_boundary():
    # American call with b == r is never exercised early -> no finite boundary.
    assert baw_critical_spot(100.0, 1.0, 0.05, 0.2, "call") is None


def test_call_boundary_above_strike_and_hits_intrinsic():
    Sc = baw_critical_spot(100.0, 1.0, 0.05, 0.2, "call", b=0.0)
    assert Sc > 100.0
    # At S* the American price equals the exercise intrinsic S* - K.
    assert baw_american(Sc, 100.0, 1.0, 0.05, 0.2, "call", b=0.0) == pytest.approx(
        Sc - 100.0, abs=1e-3)


def test_put_boundary_below_strike_and_hits_intrinsic():
    Sp = baw_critical_spot(100.0, 1.0, 0.05, 0.2, "put")
    assert Sp < 100.0
    assert baw_american(Sp, 100.0, 1.0, 0.05, 0.2, "put") == pytest.approx(
        100.0 - Sp, abs=1e-3)


@pytest.mark.parametrize("S,ot,bb", [
    (100.0, "put", None),
    (90.0, "put", None),
    (110.0, "put", None),
    (100.0, "call", 0.0),
    (110.0, "call", 0.0),
])
def test_greeks_match_finite_difference(S, ot, bb):
    g = baw_american_greeks(S, 100.0, 1.0, 0.05, 0.2, ot, b=bb)
    h = 1e-4 * S
    fd_d = (baw_american(S + h, 100.0, 1.0, 0.05, 0.2, ot, b=bb)
            - baw_american(S - h, 100.0, 1.0, 0.05, 0.2, ot, b=bb)) / (2 * h)
    assert g["delta"] == pytest.approx(fd_d, abs=1e-4)
    assert g["price"] == pytest.approx(
        baw_american(S, 100.0, 1.0, 0.05, 0.2, ot, b=bb), abs=1e-12)


def test_put_delta_and_gamma_signs():
    g = baw_american_greeks(100.0, 100.0, 1.0, 0.05, 0.2, "put")
    assert -1.0 < g["delta"] < 0.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_american_at_least_european():
    # BAW American value dominates the European price.
    assert baw_american(90.0, 100.0, 1.0, 0.05, 0.2, "put") >= put_price(
        90.0, 100.0, 1.0, 0.05, 0.2)
    assert baw_american(110.0, 100.0, 1.0, 0.05, 0.2, "call", b=0.0) >= call_price(
        110.0, 100.0, 1.0, 0.05, 0.2, b=0.0)
