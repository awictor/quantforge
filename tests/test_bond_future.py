"""Bond futures: conversion factor, basis, cheapest-to-deliver."""

import pytest

from quantforge import (
    conversion_factor, invoice_price, gross_basis, net_basis,
    implied_repo_rate, cheapest_to_deliver,
    bond_future_dv01, futures_dv01, futures_hedge_ratio,
)


def test_conversion_factor_par():
    assert conversion_factor(0.06, 10, 0.06) == pytest.approx(1.0, abs=1e-9)


def test_conversion_factor_above_below_par():
    assert conversion_factor(0.08, 10, 0.06) > 1
    assert conversion_factor(0.04, 10, 0.06) < 1


def test_invoice_price():
    assert invoice_price(120, 0.95, 1.5) == pytest.approx(120 * 0.95 + 1.5)


def test_gross_and_net_basis():
    assert gross_basis(115, 120, 0.95) == pytest.approx(115 - 120 * 0.95)
    assert net_basis(115, 120, 0.95, 0.5) == pytest.approx(115 - 120 * 0.95 - 0.5)


def test_cheapest_to_deliver_is_min_net_basis():
    bonds = [{"price": 115, "cf": 0.95, "carry": 0.5},
             {"price": 118, "cf": 0.98, "carry": 0.4},
             {"price": 112, "cf": 0.92, "carry": 0.6}]
    i = cheapest_to_deliver(bonds, 120)
    nbs = [net_basis(b["price"], 120, b["cf"], b["carry"]) for b in bonds]
    assert i == nbs.index(min(nbs))


def test_implied_repo_rises_with_futures():
    rr = implied_repo_rate(115, 1.0, 120, 0.95, 1.2, 60)
    assert implied_repo_rate(115, 1.0, 122, 0.95, 1.2, 60) > rr


def test_bond_dv01_formula():
    assert bond_future_dv01(120, 8.5) == pytest.approx(8.5 * 120 * 1e-4)


def test_futures_dv01_gears_by_conversion_factor():
    ctd = bond_future_dv01(115, 9.0)
    f = futures_dv01(ctd, 0.95)
    assert f == pytest.approx(ctd / 0.95)
    assert f > ctd   # CF < 1 gears the futures DV01 up


def test_hedge_ratio_offsets_bond_dv01():
    ctd = bond_future_dv01(115, 9.0)
    f = futures_dv01(ctd, 0.95)
    bond = bond_future_dv01(120, 8.5)
    n = futures_hedge_ratio(bond, f)
    assert n == pytest.approx(bond / f)
    assert n * f == pytest.approx(bond)


def test_hedge_ratio_scales_with_bond_dv01():
    ctd = bond_future_dv01(115, 9.0)
    f = futures_dv01(ctd, 0.95)
    bond = bond_future_dv01(120, 8.5)
    assert futures_hedge_ratio(2 * bond, f) > futures_hedge_ratio(bond, f)


def test_dv01_validation():
    with pytest.raises(ValueError):
        futures_dv01(1.0, 0)
    with pytest.raises(ValueError):
        futures_hedge_ratio(1.0, 0)


def test_validation():
    with pytest.raises(ValueError):
        conversion_factor(0.06, 0, 0.06)
    with pytest.raises(ValueError):
        cheapest_to_deliver([], 120)
    with pytest.raises(ValueError):
        implied_repo_rate(115, 1.0, 120, 0.95, 1.2, 0)
