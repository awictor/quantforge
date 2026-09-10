"""Greeks of a rough-Bergomi call by CRN bumps (rbergomi_greeks_cv)."""

import pytest

from quantforge import rbergomi_greeks_cv, rbergomi_price_cv


S, K, T = 100.0, 100.0, 1.0
XI0, ETA, H, RHO, R = 0.04, 1.5, 0.1, -0.7, 0.02


@pytest.mark.slow
def test_delta_matches_crn_finite_difference():
    g = rbergomi_greeks_cv(S, K, T, XI0, ETA, H, RHO, R, n_steps=50,
                           n_paths=20_000, seed=1)
    h = 1.0
    up = rbergomi_price_cv(S + h, K, T, XI0, ETA, H, RHO, R, 50, 20_000,
                           seed=1).price
    dn = rbergomi_price_cv(S - h, K, T, XI0, ETA, H, RHO, R, 50, 20_000,
                           seed=1).price
    fd = (up - dn) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-6)


def test_call_greek_signs():
    g = rbergomi_greeks_cv(S, K, T, XI0, ETA, H, RHO, R, n_steps=20,
                           n_paths=4_000, seed=2)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega_xi0"] > 0.0


def test_reproducible():
    kw = dict(n_steps=30, n_paths=6_000, seed=99)
    a = rbergomi_greeks_cv(S, K, T, XI0, ETA, H, RHO, R, **kw)
    b = rbergomi_greeks_cv(S, K, T, XI0, ETA, H, RHO, R, **kw)
    assert a["delta"] == b["delta"]
    assert a["price"] == b["price"]


def test_price_field_matches_price_cv():
    g = rbergomi_greeks_cv(S, K, T, XI0, ETA, H, RHO, R, n_steps=20,
                           n_paths=4_000, seed=3)
    p = rbergomi_price_cv(S, K, T, XI0, ETA, H, RHO, R, 20, 4_000, seed=3).price
    assert g["price"] == pytest.approx(p, abs=1e-9)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        rbergomi_greeks_cv(-1, K, T, XI0, ETA, H, RHO, R)
