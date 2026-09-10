"""Tests for the simple chooser option."""

import math

import pytest

from quantforge import chooser_option, call_price, put_price, OptionType


def test_choose_at_expiry_is_straddle():
    # t_choose = T: choosing at expiry is exactly a straddle.
    S, K, T, r, sigma = 100, 100, 1.0, 0.05, 0.25
    straddle = call_price(S, K, T, r, sigma) + put_price(S, K, T, r, sigma)
    assert chooser_option(S, K, T, T, r, sigma) == pytest.approx(straddle, abs=1e-9)


def test_chooser_between_vanilla_and_straddle():
    S, K, T, r, sigma = 100, 100, 1.0, 0.05, 0.25
    c = chooser_option(S, K, 0.5, T, r, sigma)
    call = call_price(S, K, T, r, sigma)
    put = put_price(S, K, T, r, sigma)
    assert c > max(call, put)          # optionality is worth more than either
    assert c < call + put              # but less than holding both (a straddle)


def test_earlier_choice_is_cheaper():
    # Choosing sooner gives up optionality, so the value is lower.
    S, K, T, r, sigma = 100, 100, 1.0, 0.05, 0.25
    early = chooser_option(S, K, 0.1, T, r, sigma)
    late = chooser_option(S, K, 0.9, T, r, sigma)
    assert late > early


def test_matches_reference_value():
    # Cross-checked against a choice-date Monte Carlo (~17.08).
    v = chooser_option(100, 100, 0.5, 1.0, 0.05, 0.25)
    assert v == pytest.approx(17.04, abs=0.1)


def test_positive_across_vols():
    for sigma in (0.1, 0.3, 0.6):
        assert chooser_option(100, 100, 0.5, 1.0, 0.05, sigma) > 0


def test_rejects_bad_choice_time():
    with pytest.raises(ValueError):
        chooser_option(100, 100, 1.5, 1.0, 0.05, 0.25)   # t_choose > T
    with pytest.raises(ValueError):
        chooser_option(100, 100, -0.1, 1.0, 0.05, 0.25)  # negative
