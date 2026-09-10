"""Greeks of a composite (compo) FX option (compo_option_greeks)."""

import pytest

from quantforge import (
    compo_option_greeks,
    compo_option,
    delta as bsm_delta,
    OptionType,
)


S, K, T = 100.0, 100.0, 1.0
RD, RF = 0.04, 0.06
SA, SFX, RHO, Q = 0.2, 0.1, 0.3, 0.02


def test_delta_matches_finite_difference():
    g = compo_option_greeks(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL)
    h = 0.01
    fd = (compo_option(S + h, K, T, RD, RF, SA, SFX, RHO, Q)
          - compo_option(S - h, K, T, RD, RF, SA, SFX, RHO, Q)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-5)


def test_zero_fx_vol_matches_plain_bsm():
    # sigma_fx = 0: the compo vol collapses to sigma_asset, carry b = rd - q.
    g = compo_option_greeks(S, K, T, RD, RF, SA, 0.0, RHO, Q, OptionType.CALL)
    bd = bsm_delta(S, K, T, RD, SA, OptionType.CALL, b=RD - Q)
    assert g["delta"] == pytest.approx(bd, abs=1e-9)


def test_compo_is_long_fx_vol():
    # Unlike a quanto, a compo option's payoff scales with the FX rate, so it is
    # long FX volatility: positive fx_vega.
    g = compo_option_greeks(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL)
    assert g["fx_vega"] > 0.0


def test_call_gamma_and_vega_positive():
    g = compo_option_greeks(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL)
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_price_field_matches_compo_option():
    g = compo_option_greeks(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL)
    assert g["price"] == pytest.approx(
        compo_option(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL),
        abs=1e-12)


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        compo_option_greeks(S, K, T, RD, RF, SA, SFX, 1.5, Q)
