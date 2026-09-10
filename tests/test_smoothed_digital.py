"""Smoothed (call-spread) pathwise delta of a cash-or-nothing digital."""

import pytest

from quantforge import (
    OptionType,
    smoothed_digital_delta,
    cash_or_nothing,
)


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


def _analytic_digital_delta():
    h = 0.01
    return (cash_or_nothing(S + h, K, T, R, SIG, OptionType.CALL)
            - cash_or_nothing(S - h, K, T, R, SIG, OptionType.CALL)) / (2 * h)


@pytest.mark.slow
def test_matches_analytic_digital_delta():
    sd = smoothed_digital_delta(S, K, T, R, SIG, OptionType.CALL,
                                eps_rel=0.02, n_paths=400_000, seed=1)
    assert sd.price == pytest.approx(_analytic_digital_delta(),
                                     abs=3.0 * sd.std_error + 1e-4)


@pytest.mark.slow
def test_narrower_spread_lowers_bias():
    an = _analytic_digital_delta()
    wide = smoothed_digital_delta(S, K, T, R, SIG, eps_rel=0.1,
                                  n_paths=400_000, seed=1)
    narrow = smoothed_digital_delta(S, K, T, R, SIG, eps_rel=0.02,
                                    n_paths=400_000, seed=1)
    assert abs(narrow.price - an) < abs(wide.price - an)


def test_positive_for_call():
    sd = smoothed_digital_delta(S, K, T, R, SIG, OptionType.CALL,
                                n_paths=100_000, seed=2)
    assert sd.price > 0.0


def test_put_digital_delta_negative():
    sd = smoothed_digital_delta(S, K, T, R, SIG, OptionType.PUT,
                                n_paths=100_000, seed=3)
    assert sd.price < 0.0


def test_bad_params_raise():
    with pytest.raises(ValueError):
        smoothed_digital_delta(-1, K, T, R, SIG)
