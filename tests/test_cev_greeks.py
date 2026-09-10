"""Greeks of a CEV option (cev_greeks)."""

import pytest

from quantforge import cev_greeks, cev_price, OptionType


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


@pytest.mark.parametrize("beta", [0.3, 0.5, 0.7])
def test_delta_matches_finite_difference(beta):
    g = cev_greeks(S, K, T, R, SIG, beta, OptionType.CALL)
    h = 0.01
    fd = (cev_price(S + h, K, T, R, SIG, beta)
          - cev_price(S - h, K, T, R, SIG, beta)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)


def test_call_greek_signs():
    g = cev_greeks(S, K, T, R, SIG, 0.5, OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_put_delta_negative():
    g = cev_greeks(S, K, T, R, SIG, 0.5, OptionType.PUT)
    assert g["delta"] < 0.0


def test_lower_beta_lifts_atm_delta():
    # A lower beta steepens the local-vol skew; the ATM call delta rises.
    g_lo = cev_greeks(S, K, T, R, SIG, 0.3, OptionType.CALL)
    g_hi = cev_greeks(S, K, T, R, SIG, 0.7, OptionType.CALL)
    assert g_lo["delta"] > g_hi["delta"]


def test_price_field_matches_price():
    g = cev_greeks(S, K, T, R, SIG, 0.5, OptionType.CALL)
    assert g["price"] == pytest.approx(
        cev_price(S, K, T, R, SIG, 0.5, OptionType.CALL), abs=1e-12)


def test_bad_beta_raises():
    with pytest.raises(ValueError):
        cev_greeks(S, K, T, R, SIG, 1.0)
