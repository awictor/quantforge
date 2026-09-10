"""Greeks of a forward-start option (forward_start_greeks)."""

import pytest

from quantforge import forward_start_greeks, forward_start_price, OptionType


S, TS, TE, R, SIG = 100.0, 0.5, 1.0, 0.05, 0.2


def test_delta_matches_finite_difference():
    g = forward_start_greeks(S, TS, TE, R, SIG, 1.0, OptionType.CALL)
    h = 0.01
    fd = (forward_start_price(S + h, TS, TE, R, SIG)
          - forward_start_price(S - h, TS, TE, R, SIG)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-6)


def test_gamma_is_zero():
    # The price is exactly linear in the spot before the strike is fixed.
    g = forward_start_greeks(S, TS, TE, R, SIG, 1.0, OptionType.CALL)
    assert g["gamma"] == 0.0


def test_delta_constant_in_spot():
    # Linear in S -> delta does not depend on the spot level.
    g1 = forward_start_greeks(100.0, TS, TE, R, SIG)
    g2 = forward_start_greeks(150.0, TS, TE, R, SIG)
    assert g1["delta"] == pytest.approx(g2["delta"], abs=1e-12)


def test_vega_matches_finite_difference():
    g = forward_start_greeks(S, TS, TE, R, SIG, 1.0, OptionType.CALL)
    h = 1e-4
    fd = (forward_start_price(S, TS, TE, R, SIG + h)
          - forward_start_price(S, TS, TE, R, SIG - h)) / (2 * h)
    assert g["vega"] == pytest.approx(fd, abs=1e-3)


def test_price_field_matches_price():
    g = forward_start_greeks(S, TS, TE, R, SIG, 1.0, OptionType.CALL)
    assert g["price"] == pytest.approx(
        forward_start_price(S, TS, TE, R, SIG, 1.0, OptionType.CALL), abs=1e-12)


def test_theta_nonzero_under_carry():
    # With b != r the discounted scale factor decays, so theta is nonzero.
    g = forward_start_greeks(S, TS, TE, R, SIG, 1.0, OptionType.CALL, b=0.02)
    assert abs(g["theta"]) > 1e-3


def test_bad_times_raise():
    with pytest.raises(ValueError):
        forward_start_greeks(S, 1.0, 0.5, R, SIG)  # t_start >= t_expiry
