"""Greeks of a gap option (gap_option_greeks)."""

import pytest

from quantforge import gap_option_greeks, gap_option, greeks, OptionType


S, T, R, SIG = 100.0, 1.0, 0.05, 0.2


def test_equal_strikes_reduce_to_vanilla():
    g = gap_option_greeks(S, 100.0, 100.0, T, R, SIG, OptionType.CALL)
    van = greeks(S, 100.0, T, R, SIG, OptionType.CALL)
    assert g["delta"] == pytest.approx(van.delta, abs=1e-4)
    assert g["gamma"] == pytest.approx(van.gamma, abs=1e-5)
    assert g["vega"] == pytest.approx(van.vega, abs=1e-2)


def test_delta_matches_finite_difference():
    g = gap_option_greeks(S, 110.0, 100.0, T, R, SIG, OptionType.CALL)
    h = 0.01
    fd = (gap_option(S + h, 110.0, 100.0, T, R, SIG)
          - gap_option(S - h, 110.0, 100.0, T, R, SIG)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-5)


def test_price_field_matches_gap_option():
    g = gap_option_greeks(S, 110.0, 100.0, T, R, SIG, OptionType.CALL)
    assert g["price"] == pytest.approx(
        gap_option(S, 110.0, 100.0, T, R, SIG, OptionType.CALL), abs=1e-12)


def test_put_matches_finite_difference():
    g = gap_option_greeks(S, 90.0, 100.0, T, R, SIG, OptionType.PUT)
    h = 0.01
    fd = (gap_option(S + h, 90.0, 100.0, T, R, SIG, OptionType.PUT)
          - gap_option(S - h, 90.0, 100.0, T, R, SIG, OptionType.PUT)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-5)


def test_bad_payoff_strike_raises():
    with pytest.raises(ValueError):
        gap_option_greeks(S, 110.0, -1.0, T, R, SIG)
