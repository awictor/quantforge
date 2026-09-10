"""Cash-or-nothing digital Greeks by the likelihood-ratio method."""

import pytest

from quantforge import (
    OptionType,
    lr_digital_greeks,
    digital_greeks,
    cash_or_nothing,
)


S, K, T, R, SIG = 100.0, 100.0, 0.5, 0.05, 0.2


def _analytic_vega(ot):
    h = 1e-4
    up = cash_or_nothing(S, K, T, R, SIG + h, ot)
    dn = cash_or_nothing(S, K, T, R, SIG - h, ot)
    return (up - dn) / (2 * h)


@pytest.mark.slow
def test_call_greeks_match_analytic():
    an = digital_greeks(S, K, T, R, SIG, OptionType.CALL)
    g = lr_digital_greeks(S, K, T, R, SIG, OptionType.CALL,
                          n_paths=800_000, seed=1)
    assert g["price"] == pytest.approx(an["price"], abs=3.0 * g["price_se"])
    assert g["delta"] == pytest.approx(an["delta"], abs=3.0 * g["delta_se"])
    assert g["gamma"] == pytest.approx(an["gamma"], abs=3.0 * g["gamma_se"])
    assert g["vega"] == pytest.approx(_analytic_vega(OptionType.CALL),
                                      abs=3.0 * g["vega_se"])


@pytest.mark.slow
def test_put_greeks_match_analytic():
    an = digital_greeks(S, K, T, R, SIG, OptionType.PUT)
    g = lr_digital_greeks(S, K, T, R, SIG, OptionType.PUT,
                          n_paths=800_000, seed=2)
    assert g["delta"] == pytest.approx(an["delta"], abs=3.0 * g["delta_se"])
    assert g["gamma"] == pytest.approx(an["gamma"], abs=3.0 * g["gamma_se"])
    assert g["vega"] == pytest.approx(_analytic_vega(OptionType.PUT),
                                      abs=3.0 * g["vega_se"])


def test_call_delta_positive_put_delta_negative():
    gc = lr_digital_greeks(S, K, T, R, SIG, OptionType.CALL,
                           n_paths=100_000, seed=3)
    gp = lr_digital_greeks(S, K, T, R, SIG, OptionType.PUT,
                           n_paths=100_000, seed=4)
    assert gc["delta"] > 0.0
    assert gp["delta"] < 0.0


def test_price_matches_closed_form():
    g = lr_digital_greeks(S, K, T, R, SIG, OptionType.CALL,
                          n_paths=200_000, seed=5)
    cf = cash_or_nothing(S, K, T, R, SIG, OptionType.CALL)
    assert g["price"] == pytest.approx(cf, abs=3.0 * g["price_se"])


def test_bad_params_raise():
    with pytest.raises(ValueError):
        lr_digital_greeks(-1, K, T, R, SIG)
