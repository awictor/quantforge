"""Greeks of a rough-Heston option (rough_heston_greeks)."""

import pytest

from quantforge import (
    rough_heston_greeks,
    rough_heston_price,
    heston_price,
    OptionType,
)


S, K, T, R = 100.0, 100.0, 1.0, 0.03
V0, KAPPA, THETA, NU, RHO = 0.04, 1.5, 0.04, 0.4, -0.6


@pytest.mark.slow
def test_h_half_matches_classical_heston():
    # H = 0.5 reduces rough-Heston to classical Heston with xi = kappa * nu.
    g = rough_heston_greeks(S, K, T, R, V0, KAPPA, THETA, NU, RHO, H=0.5)
    h = 0.5
    hd = (heston_price(S + h, K, T, R, V0, KAPPA, THETA, KAPPA * NU, RHO,
                       OptionType.CALL)
          - heston_price(S - h, K, T, R, V0, KAPPA, THETA, KAPPA * NU, RHO,
                         OptionType.CALL)) / (2 * h)
    assert g["delta"] == pytest.approx(hd, abs=2e-3)


@pytest.mark.slow
def test_rough_delta_matches_finite_difference():
    g = rough_heston_greeks(S, K, T, R, V0, KAPPA, THETA, NU, RHO, H=0.3)
    h = 0.5
    fd = (rough_heston_price(S + h, K, T, R, V0, KAPPA, THETA, NU, RHO, H=0.3,
                             n_grid=200)
          - rough_heston_price(S - h, K, T, R, V0, KAPPA, THETA, NU, RHO, H=0.3,
                               n_grid=200)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=2e-3)


@pytest.mark.slow
def test_rough_greek_signs():
    g = rough_heston_greeks(S, K, T, R, V0, KAPPA, THETA, NU, RHO, H=0.3)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega_v0"] > 0.0


def test_bad_hurst_raises():
    with pytest.raises(ValueError):
        rough_heston_greeks(S, K, T, R, V0, KAPPA, THETA, NU, RHO, H=0.7)
