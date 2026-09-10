"""Greeks of a cliquet (ratchet) option (cliquet_greeks)."""

import pytest

from quantforge import cliquet_greeks, cliquet_price, OptionType


S, R, SIG = 100.0, 0.05, 0.2
RT = [0.25, 0.5, 0.75, 1.0]


def test_delta_matches_finite_difference():
    g = cliquet_greeks(S, RT, R, SIG, 1.0, OptionType.CALL)
    h = 0.01
    fd = (cliquet_price(S + h, RT, R, SIG) - cliquet_price(S - h, RT, R, SIG)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-6)


def test_delta_equals_price_over_spot():
    # Every strike scales with S (alpha * S_reset), so the cliquet is homogeneous
    # of degree 1 in the spot: delta = price / S exactly.
    g = cliquet_greeks(S, RT, R, SIG, 1.0, OptionType.CALL)
    p = cliquet_price(S, RT, R, SIG, 1.0, OptionType.CALL)
    assert g["delta"] == pytest.approx(p / S, abs=1e-5)


def test_gamma_is_zero():
    # Homogeneous degree 1 -> no spot convexity.
    g = cliquet_greeks(S, RT, R, SIG, 1.0, OptionType.CALL)
    assert g["gamma"] == pytest.approx(0.0, abs=1e-6)


def test_vega_matches_finite_difference():
    g = cliquet_greeks(S, RT, R, SIG, 1.0, OptionType.CALL)
    h = 1e-4
    fd = (cliquet_price(S, RT, R, SIG + h) - cliquet_price(S, RT, R, SIG - h)) / (2 * h)
    assert g["vega"] == pytest.approx(fd, abs=1e-2)


def test_price_field_matches_price():
    g = cliquet_greeks(S, RT, R, SIG, 1.0, OptionType.CALL)
    assert g["price"] == pytest.approx(
        cliquet_price(S, RT, R, SIG, 1.0, OptionType.CALL), abs=1e-12)


def test_bad_reset_times_raise():
    with pytest.raises(ValueError):
        cliquet_greeks(S, [0.5, 0.5, 1.0], R, SIG)
