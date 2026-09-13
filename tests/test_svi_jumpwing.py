"""SVI raw <-> jump-wing parameterization."""

import math

import pytest

from quantforge import SVIParams, raw_to_jumpwing, jumpwing_to_raw


def _max_w_diff(p1, p2, ks):
    return max(abs(p1.total_variance(k) - p2.total_variance(k)) for k in ks)


KS = [-0.6, -0.3, -0.1, 0.0, 0.1, 0.3, 0.6]


def test_jw_v_is_atm_variance():
    p = SVIParams(a=0.04, b=0.4, rho=-0.4, m=0.05, s=0.2)
    t = 1.5
    jw = raw_to_jumpwing(p, t)
    assert abs(jw.v * t - p.total_variance(0.0)) < 1e-12


def test_jw_psi_matches_finite_difference_skew():
    p = SVIParams(a=0.04, b=0.4, rho=-0.4, m=0.05, s=0.2)
    jw = raw_to_jumpwing(p, 1.0)
    h = 1e-6
    dwdk = (p.total_variance(h) - p.total_variance(-h)) / (2 * h)
    psi_fd = dwdk / (2 * math.sqrt(p.total_variance(0.0)))
    assert abs(jw.psi - psi_fd) < 1e-6


def test_jw_wing_slopes_match_asymptotics():
    p = SVIParams(a=0.04, b=0.4, rho=-0.4, m=0.05, s=0.2)
    jw = raw_to_jumpwing(p, 1.0)
    sqrt_w = math.sqrt(p.total_variance(0.0))
    assert abs(jw.c - p.b * (1 + p.rho) / sqrt_w) < 1e-12
    assert abs(jw.p - p.b * (1 - p.rho) / sqrt_w) < 1e-12


def test_roundtrip_asymmetric():
    p = SVIParams(a=0.04, b=0.4, rho=-0.4, m=0.05, s=0.2)
    p2 = jumpwing_to_raw(raw_to_jumpwing(p, 1.0), 1.0)
    for name in ("a", "b", "rho", "m", "s"):
        assert abs(getattr(p, name) - getattr(p2, name)) < 1e-9
    assert _max_w_diff(p, p2, KS) < 1e-12


def test_roundtrip_positive_skew():
    p = SVIParams(a=0.02, b=0.3, rho=0.35, m=-0.04, s=0.18)
    p2 = jumpwing_to_raw(raw_to_jumpwing(p, 2.0), 2.0)
    assert _max_w_diff(p, p2, KS) < 1e-10


def test_roundtrip_zero_shift_nonzero_skew():
    # m = 0 with rho != 0 is recoverable.
    p = SVIParams(a=0.03, b=0.3, rho=-0.3, m=0.0, s=0.15)
    p2 = jumpwing_to_raw(raw_to_jumpwing(p, 1.0), 1.0)
    assert _max_w_diff(p, p2, KS) < 1e-9


def test_roundtrip_various_expiries():
    p = SVIParams(a=0.05, b=0.5, rho=-0.6, m=0.1, s=0.25)
    for t in (0.25, 0.5, 1.0, 3.0):
        p2 = jumpwing_to_raw(raw_to_jumpwing(p, t), t)
        assert _max_w_diff(p, p2, KS) < 1e-9


def test_validation():
    p = SVIParams(a=0.04, b=0.4, rho=-0.4, m=0.05, s=0.2)
    with pytest.raises(ValueError):
        raw_to_jumpwing(p, 0.0)
    jw = raw_to_jumpwing(p, 1.0)
    with pytest.raises(ValueError):
        jumpwing_to_raw(jw, -1.0)
