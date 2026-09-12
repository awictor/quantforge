"""AR order selection by AIC / BIC."""

import random

import pytest

from quantforge import select_ar_order, ar_information_criteria


def _ar1(n, seed):
    rng = random.Random(seed)
    y = [0.0]
    for _ in range(n):
        y.append(0.6 * y[-1] + rng.gauss(0, 1))
    return y[1:]


def _ar2(n, seed):
    rng = random.Random(seed)
    y = [0.0, 0.0]
    for _ in range(n):
        y.append(0.5 * y[-1] - 0.3 * y[-2] + rng.gauss(0, 1))
    return y[2:]


def test_bic_selects_true_order_ar1():
    order, _ = select_ar_order(_ar1(5000, 1), 8, "bic")
    assert order == 1


def test_bic_selects_true_order_ar2():
    order, _ = select_ar_order(_ar2(5000, 2), 8, "bic")
    assert order == 2


def test_bic_order_not_larger_than_aic():
    y = _ar2(5000, 2)
    o_bic, _ = select_ar_order(y, 8, "bic")
    o_aic, _ = select_ar_order(y, 8, "aic")
    assert o_bic <= o_aic


def test_aic_drops_at_true_order():
    y = _ar2(5000, 2)
    aic1 = ar_information_criteria(y, 1)[0]
    aic2 = ar_information_criteria(y, 2)[0]
    assert aic2 < aic1


def test_validation():
    y = _ar1(200, 1)
    with pytest.raises(ValueError):
        select_ar_order(y, 0)
    with pytest.raises(ValueError):
        select_ar_order([1.0, 2.0], 5)
    with pytest.raises(ValueError):
        select_ar_order(y, 5, "xic")
