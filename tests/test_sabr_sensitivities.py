"""Exact SABR parameter sensitivities (dual-number AD) vs finite differences."""

import pytest

from quantforge import (
    sabr_vol,
    sabr_sensitivities,
    sabr_jacobian,
)


F, T, BETA = 100.0, 1.0, 0.5
ALPHA, RHO, NU = 2.0, -0.3, 0.4


def _fd(var, K, h):
    kw = dict(F=F, K=K, t=T, alpha=ALPHA, beta=BETA, rho=RHO, nu=NU)
    up, dn = dict(kw), dict(kw)
    up[var] += h
    dn[var] -= h
    return (sabr_vol(**up) - sabr_vol(**dn)) / (2.0 * h)


@pytest.mark.parametrize("K", [80.0, 90.0, 110.0, 125.0])
def test_ad_partials_match_finite_difference(K):
    s = sabr_sensitivities(F, K, T, ALPHA, BETA, RHO, NU)
    assert s["vol"] == pytest.approx(sabr_vol(F, K, T, ALPHA, BETA, RHO, NU),
                                     abs=1e-14)
    checks = [
        ("F", "d_dF", 1e-3),
        ("K", "d_dK", 1e-3),
        ("alpha", "d_dalpha", 1e-6),
        ("rho", "d_drho", 1e-6),
        ("nu", "d_dnu", 1e-6),
    ]
    for var, key, h in checks:
        assert s[key] == pytest.approx(_fd(var, K, h), abs=1e-6)


def test_atm_alpha_partial_is_vol_over_alpha():
    # At the money the leading term is alpha / F^{1-beta}, so d vol / d alpha is
    # ~ vol / alpha (the higher-order time terms shift it only slightly).
    s = sabr_sensitivities(F, F, T, ALPHA, BETA, RHO, NU)
    assert s["d_dalpha"] == pytest.approx(s["vol"] / ALPHA, rel=0.05)


def test_jacobian_shape_and_columns_match_sensitivities():
    strikes = [80.0, 100.0, 120.0]
    J = sabr_jacobian(F, T, strikes, ALPHA, BETA, RHO, NU)
    assert len(J) == 3 and all(len(row) == 3 for row in J)
    for row, K in zip(J, strikes):
        s = sabr_sensitivities(F, K, T, ALPHA, BETA, RHO, NU)
        assert row == pytest.approx([s["d_dalpha"], s["d_drho"], s["d_dnu"]])


def test_vol_of_vol_raises_smile_wings():
    # d vol / d nu should be positive in the wings (more vol-of-vol = fatter
    # smile away from the money).
    lo = sabr_sensitivities(F, 80.0, T, ALPHA, BETA, RHO, NU)["d_dnu"]
    hi = sabr_sensitivities(F, 125.0, T, ALPHA, BETA, RHO, NU)["d_dnu"]
    assert lo > 0 and hi > 0


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        sabr_sensitivities(-1, 100, T, ALPHA, BETA, RHO, NU)
    with pytest.raises(ValueError):
        sabr_sensitivities(F, 100, T, 0.0, BETA, RHO, NU)
    with pytest.raises(ValueError):
        sabr_sensitivities(F, 100, 0.0, ALPHA, BETA, RHO, NU)
