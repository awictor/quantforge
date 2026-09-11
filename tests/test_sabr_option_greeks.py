"""Greeks of an option priced at the SABR smile vol (sabr_option_greeks)."""

import pytest

from quantforge import sabr_option_greeks, sabr_vol, OptionType
from quantforge.bsm import price as bp


F, K, T = 100.0, 105.0, 1.0
ALPHA, BETA, RHO, NU = 0.2, 0.5, -0.3, 0.4


def _sabr_price(f, ot):
    return bp(f, K, T, 0.0, sabr_vol(f, K, T, ALPHA, BETA, RHO, NU), ot, b=0.0)


def test_total_delta_matches_finite_difference():
    # The total delta includes the smile backbone dsigma/dF, so it must match a
    # finite difference that re-computes the SABR vol at each bumped forward.
    g = sabr_option_greeks(F, K, T, ALPHA, BETA, RHO, NU, OptionType.CALL)
    h = 0.01
    fd = (_sabr_price(F + h, OptionType.CALL)
          - _sabr_price(F - h, OptionType.CALL)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)


def test_total_delta_differs_from_black_delta():
    # The backbone term makes the SABR delta differ from the vol-fixed Black
    # delta (nonzero for beta != 1 / rho != 0).
    g = sabr_option_greeks(F, K, T, ALPHA, BETA, RHO, NU, OptionType.CALL)
    assert abs(g["delta"] - g["black_delta"]) > 1e-4


def test_vega_positive():
    g = sabr_option_greeks(F, K, T, ALPHA, BETA, RHO, NU, OptionType.CALL)
    assert g["vega"] > 0.0


def test_vol_matches_sabr_vol():
    g = sabr_option_greeks(F, K, T, ALPHA, BETA, RHO, NU, OptionType.CALL)
    assert g["vol"] == pytest.approx(sabr_vol(F, K, T, ALPHA, BETA, RHO, NU),
                                     abs=1e-12)


def test_discount_scales_price_and_greeks():
    g1 = sabr_option_greeks(F, K, T, ALPHA, BETA, RHO, NU, OptionType.CALL,
                            discount=1.0)
    g2 = sabr_option_greeks(F, K, T, ALPHA, BETA, RHO, NU, OptionType.CALL,
                            discount=0.95)
    assert g2["price"] == pytest.approx(0.95 * g1["price"], rel=1e-9)
    assert g2["delta"] == pytest.approx(0.95 * g1["delta"], rel=1e-9)


def test_put_price_field_matches_black():
    g = sabr_option_greeks(F, K, T, ALPHA, BETA, RHO, NU, OptionType.PUT)
    assert g["price"] == pytest.approx(_sabr_price(F, OptionType.PUT), abs=1e-12)
