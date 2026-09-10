"""Tests for the theta/carry decomposition report."""

import pytest

from quantforge import Contract, theta_carry_report, price_book


def _leg(qty=1, b=0.0, r=0.0):
    return Contract(S=100, K=100, t=0.5, r=r, sigma=0.25, option_type="call",
                    qty=qty, b=b)


def test_zero_carry_theta_equals_gamma_rent():
    # b = 0, r = 0: time decay is pure gamma rent, residual ~ 0.
    rep = theta_carry_report([_leg(b=0.0, r=0.0)])
    assert rep.residual == pytest.approx(0.0, abs=1e-6)
    assert rep.theta == pytest.approx(rep.gamma_rent, abs=1e-6)


def test_decomposition_adds_up():
    rep = theta_carry_report([_leg(b=0.05, r=0.05)])
    assert rep.theta == pytest.approx(rep.gamma_rent + rep.residual, abs=1e-9)


def test_financing_residual_nonzero_with_rate():
    rep = theta_carry_report([_leg(b=0.05, r=0.05)])
    assert rep.residual != pytest.approx(0.0, abs=1e-3)


def test_long_gamma_pays_rent_short_collects():
    long = theta_carry_report([_leg(qty=1, b=0.0, r=0.0)])
    short = theta_carry_report([_leg(qty=-1, b=0.0, r=0.0)])
    assert long.theta < 0        # long options bleed theta
    assert short.theta > 0       # short options collect it
    assert short.gamma_rent == pytest.approx(-long.gamma_rent)


def test_net_theta_matches_price_book():
    legs = [_leg(qty=3, b=0.04, r=0.04), _leg(qty=-1, b=0.04, r=0.04)]
    rep = theta_carry_report(legs)
    net = price_book(legs).net.theta
    assert rep.theta == pytest.approx(net, abs=1e-9)


def test_scales_with_quantity():
    one = theta_carry_report([_leg(qty=1)])
    ten = theta_carry_report([_leg(qty=10)])
    assert ten.gamma_rent == pytest.approx(10 * one.gamma_rent)
