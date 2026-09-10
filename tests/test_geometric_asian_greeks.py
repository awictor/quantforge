"""Greeks of the continuously-monitored geometric-average Asian (Kemna-Vorst)."""

import pytest

from quantforge import geometric_asian_greeks, geometric_asian, OptionType


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


def _fd(f, x, h):
    return (f(x + h) - f(x - h)) / (2 * h)


def test_delta_matches_finite_difference():
    g = geometric_asian_greeks(S, K, T, R, SIG, OptionType.CALL)
    fd = _fd(lambda s: geometric_asian(s, K, T, R, SIG), S, 0.01)
    assert g["delta"] == pytest.approx(fd, abs=1e-6)


def test_gamma_matches_finite_difference():
    g = geometric_asian_greeks(S, K, T, R, SIG, OptionType.CALL)
    h = 0.01
    fg = (geometric_asian(S + h, K, T, R, SIG)
          - 2 * geometric_asian(S, K, T, R, SIG)
          + geometric_asian(S - h, K, T, R, SIG)) / (h * h)
    assert g["gamma"] == pytest.approx(fg, abs=1e-5)


def test_vega_matches_finite_difference():
    g = geometric_asian_greeks(S, K, T, R, SIG, OptionType.CALL)
    fd = _fd(lambda s: geometric_asian(S, K, T, R, s), SIG, 1e-4)
    assert g["vega"] == pytest.approx(fd, abs=1e-2)


def test_call_greek_signs():
    g = geometric_asian_greeks(S, K, T, R, SIG, OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_put_delta_negative():
    g = geometric_asian_greeks(S, K, T, R, SIG, OptionType.PUT)
    assert g["delta"] < 0.0


def test_price_field_matches_geometric_asian():
    g = geometric_asian_greeks(S, K, T, R, SIG, OptionType.CALL)
    assert g["price"] == pytest.approx(
        geometric_asian(S, K, T, R, SIG, OptionType.CALL), abs=1e-12)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        geometric_asian_greeks(-1, K, T, R, SIG)
