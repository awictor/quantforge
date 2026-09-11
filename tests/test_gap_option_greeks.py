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


@pytest.mark.parametrize("Kt,Kp,ot", [
    (95.0, 105.0, OptionType.CALL),
    (110.0, 90.0, OptionType.CALL),
    (95.0, 105.0, OptionType.PUT),
    (110.0, 90.0, OptionType.PUT),
])
def test_gamma_matches_finite_difference(Kt, Kp, ot):
    g = gap_option_greeks(S, Kt, Kp, T, R, SIG, ot)
    h = 1e-4 * S
    fd = (gap_option(S + h, Kt, Kp, T, R, SIG, ot)
          - 2 * gap_option(S, Kt, Kp, T, R, SIG, ot)
          + gap_option(S - h, Kt, Kp, T, R, SIG, ot)) / (h * h)
    assert g["gamma"] == pytest.approx(fd, abs=1e-4)


def test_gamma_same_for_call_and_put():
    gc = gap_option_greeks(S, 110.0, 100.0, T, R, SIG, OptionType.CALL)
    gp = gap_option_greeks(S, 110.0, 100.0, T, R, SIG, OptionType.PUT)
    assert gc["gamma"] == pytest.approx(gp["gamma"], abs=1e-12)


def test_equal_strikes_gamma_matches_vanilla_put():
    g = gap_option_greeks(S, 100.0, 100.0, T, R, SIG, OptionType.PUT)
    van = greeks(S, 100.0, T, R, SIG, OptionType.PUT)
    assert g["delta"] == pytest.approx(van.delta, abs=1e-4)
    assert g["gamma"] == pytest.approx(van.gamma, abs=1e-5)


def test_bad_payoff_strike_raises():
    with pytest.raises(ValueError):
        gap_option_greeks(S, 110.0, -1.0, T, R, SIG)
