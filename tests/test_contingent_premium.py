"""Pay-later / contingent-premium options."""

import pytest

from quantforge import (
    contingent_premium_option, pay_later_option_value, cash_or_nothing, OptionType,
)
from quantforge.bsm import price as bsm_price


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.25


def test_fair_premium_exceeds_vanilla():
    van = bsm_price(S, K, T, R, SIG, OptionType.CALL)
    assert contingent_premium_option(S, K, T, R, SIG) > van


def test_pay_later_zero_at_fair_premium():
    prem = contingent_premium_option(S, K, T, R, SIG)
    assert pay_later_option_value(S, K, T, R, SIG, prem) == pytest.approx(0.0, abs=1e-9)


def test_pay_later_sign_around_fair():
    prem = contingent_premium_option(S, K, T, R, SIG)
    assert pay_later_option_value(S, K, T, R, SIG, prem * 0.5) > 0
    assert pay_later_option_value(S, K, T, R, SIG, prem * 1.5) < 0


def test_premium_formula():
    van = bsm_price(S, K, T, R, SIG, OptionType.CALL)
    unit = cash_or_nothing(S, K, T, R, SIG, OptionType.CALL, cash=1.0)
    assert contingent_premium_option(S, K, T, R, SIG) == pytest.approx(van / unit)


def test_deep_itm_premium_approaches_vanilla():
    deep = contingent_premium_option(200, 100, T, R, SIG) / bsm_price(
        200, 100, T, R, SIG, OptionType.CALL)
    atm = contingent_premium_option(S, K, T, R, SIG) / bsm_price(
        S, K, T, R, SIG, OptionType.CALL)
    assert abs(deep - 1) < abs(atm - 1)


def test_put_fair_premium_exceeds_vanilla():
    assert contingent_premium_option(S, K, T, R, SIG, OptionType.PUT) > \
        bsm_price(S, K, T, R, SIG, OptionType.PUT)
