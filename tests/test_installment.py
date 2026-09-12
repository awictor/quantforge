"""Installment option on a binomial tree."""

import pytest

from quantforge import installment_call, american_price

S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.25


def _tree_vanilla():
    return american_price(S, K, T, R, SIG, option_type="call", b=R, steps=200)


def test_zero_installment_equals_tree_vanilla():
    assert abs(installment_call(S, K, T, R, SIG, 0.0, [0.5], 200) - _tree_vanilla()) < 1e-9


def test_empty_paytimes_equals_vanilla():
    assert abs(installment_call(S, K, T, R, SIG, 5.0, [], 200) - _tree_vanilla()) < 1e-9


def test_higher_installment_lowers_upfront():
    lo = installment_call(S, K, T, R, SIG, 2.0, [0.5], 200)
    hi = installment_call(S, K, T, R, SIG, 5.0, [0.5], 200)
    assert hi < lo


def test_paying_costs_less_than_vanilla():
    assert installment_call(S, K, T, R, SIG, 2.0, [0.5], 200) < _tree_vanilla()


def test_prohibitive_installment_drives_to_zero():
    assert installment_call(S, K, T, R, SIG, 1e6, [0.5], 200) == 0.0


def test_more_payment_dates_lower_value():
    one = installment_call(S, K, T, R, SIG, 2.0, [0.5], 200)
    three = installment_call(S, K, T, R, SIG, 2.0, [0.25, 0.5, 0.75], 200)
    assert three < one


def test_value_nonnegative():
    assert installment_call(S, K, T, R, SIG, 10.0, [0.25, 0.5, 0.75], 200) >= 0.0


def test_validation():
    with pytest.raises(ValueError):
        installment_call(-1, K, T, R, SIG, 2.0, [0.5])
    with pytest.raises(ValueError):
        installment_call(S, K, T, R, SIG, -1.0, [0.5])
    with pytest.raises(ValueError):
        installment_call(S, K, T, R, SIG, 2.0, [1.5])  # outside (0, t)
    with pytest.raises(ValueError):
        installment_call(S, K, T, R, SIG, 2.0, [0.5, 0.5])  # duplicate
