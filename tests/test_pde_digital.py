"""Digital and no-touch binaries via the Crank-Nicolson PDE."""

import math

import pytest

from quantforge import (
    OptionType,
    crank_nicolson_digital,
    crank_nicolson_no_touch,
    cash_or_nothing,
    no_touch,
    one_touch,
)


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.25


@pytest.mark.slow
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_digital_matches_closed_form(ot):
    pde = crank_nicolson_digital(S, K, T, R, SIGMA, ot,
                                 n_space=600, n_time=600)
    cf = cash_or_nothing(S, K, T, R, SIGMA, ot)
    assert pde == pytest.approx(cf, abs=1e-2)


def test_digital_call_plus_put_is_discount():
    c = crank_nicolson_digital(S, K, T, R, SIGMA, OptionType.CALL,
                               n_space=500, n_time=500)
    p = crank_nicolson_digital(S, K, T, R, SIGMA, OptionType.PUT,
                               n_space=500, n_time=500)
    assert (c + p) == pytest.approx(math.exp(-R * T), abs=1e-3)


@pytest.mark.slow
@pytest.mark.parametrize("H", [85, 120])
def test_no_touch_matches_closed_form(H):
    pde = crank_nicolson_no_touch(S, H, T, R, SIGMA, n_space=800, n_time=800)
    cf = no_touch(S, H, T, R, SIGMA)
    assert pde == pytest.approx(cf, abs=1.5e-2)


@pytest.mark.slow
def test_one_touch_complement():
    # pay-at-expiry one-touch = disc*cash - no_touch.
    nt = crank_nicolson_no_touch(S, 85, T, R, SIGMA, n_space=800, n_time=800)
    ot_pay = math.exp(-R * T) - nt
    assert ot_pay == pytest.approx(one_touch(S, 85, T, R, SIGMA,
                                             payoff_at_hit=False), abs=1.5e-2)


def test_no_touch_below_cash():
    nt = crank_nicolson_no_touch(S, 85, T, R, SIGMA, n_space=400, n_time=400)
    assert 0.0 < nt < 1.0


def test_bad_params_raise():
    with pytest.raises(ValueError):
        crank_nicolson_digital(S, K, T, R)                # no vol
    with pytest.raises(ValueError):
        crank_nicolson_no_touch(S, -1, T, R, SIGMA)       # bad barrier
