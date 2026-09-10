"""Greeks of a simple chooser option (chooser_option_greeks)."""

import pytest

from quantforge import chooser_option_greeks, chooser_option


S, K, TC, T, R, SIG = 100.0, 100.0, 0.5, 1.0, 0.05, 0.2


def _fd(f, x, h):
    return (f(x + h) - f(x - h)) / (2 * h)


def test_delta_matches_finite_difference():
    g = chooser_option_greeks(S, K, TC, T, R, SIG)
    fd = _fd(lambda s: chooser_option(s, K, TC, T, R, SIG), S, 0.01)
    assert g["delta"] == pytest.approx(fd, abs=1e-6)


def test_gamma_matches_finite_difference():
    g = chooser_option_greeks(S, K, TC, T, R, SIG)
    h = 0.01
    fg = (chooser_option(S + h, K, TC, T, R, SIG)
          - 2 * chooser_option(S, K, TC, T, R, SIG)
          + chooser_option(S - h, K, TC, T, R, SIG)) / (h * h)
    assert g["gamma"] == pytest.approx(fg, abs=1e-5)


def test_vega_matches_finite_difference():
    g = chooser_option_greeks(S, K, TC, T, R, SIG)
    fd = _fd(lambda s: chooser_option(S, K, TC, T, R, s), SIG, 1e-4)
    assert g["vega"] == pytest.approx(fd, abs=1e-3)


def test_gamma_and_vega_positive():
    # Long both a call and a put -> convex, positive gamma and vega.
    g = chooser_option_greeks(S, K, TC, T, R, SIG)
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_price_field_matches_chooser_option():
    g = chooser_option_greeks(S, K, TC, T, R, SIG)
    assert g["price"] == pytest.approx(chooser_option(S, K, TC, T, R, SIG),
                                       abs=1e-12)


def test_bad_choose_time_raises():
    with pytest.raises(ValueError):
        chooser_option_greeks(S, K, 1.5, T, R, SIG)
