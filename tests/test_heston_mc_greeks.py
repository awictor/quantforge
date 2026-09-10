"""Heston Greeks by common-random-number bumps on the QE MC (heston_mc_greeks)."""

import pytest

from quantforge import heston_mc_greeks, heston_price, OptionType


S, K, T, R = 100.0, 100.0, 1.0, 0.03
V0, KAPPA, THETA, XI, RHO = 0.04, 1.5, 0.04, 0.5, -0.7


def _fd(f, x, h):
    return (f(x + h) - f(x - h)) / (2 * h)


def _P(**over):
    kw = dict(S=S, v0=V0, theta=THETA, xi=XI, rho=RHO)
    kw.update(over)
    return heston_price(kw["S"], K, T, R, kw["v0"], KAPPA, kw["theta"],
                        kw["xi"], kw["rho"], OptionType.CALL)


@pytest.mark.slow
def test_delta_gamma_match_fourier():
    g = heston_mc_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL,
                         n_steps=60, n_paths=120_000, seed=1)
    assert g["delta"] == pytest.approx(_fd(lambda s: _P(S=s), S, 0.5), abs=5e-3)
    d2 = (_P(S=S + 0.5) - 2 * _P() + _P(S=S - 0.5)) / 0.25
    assert g["gamma"] == pytest.approx(d2, abs=3e-3)


@pytest.mark.slow
def test_variance_vegas_match_fourier():
    g = heston_mc_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL,
                         n_steps=60, n_paths=120_000, seed=2)
    assert g["vega_v0"] == pytest.approx(
        _fd(lambda v: _P(v0=v), V0, 1e-3), rel=0.02)
    assert g["vega_theta"] == pytest.approx(
        _fd(lambda x: _P(theta=x), THETA, 1e-3), rel=0.02)


@pytest.mark.slow
def test_volvol_and_rho_match_fourier():
    g = heston_mc_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL,
                         n_steps=60, n_paths=180_000, seed=3)
    assert g["volvol"] == pytest.approx(_fd(lambda x: _P(xi=x), XI, 1e-3), abs=0.3)
    assert g["rho_sens"] == pytest.approx(
        _fd(lambda x: _P(rho=x), RHO, 1e-3), abs=0.1)


def test_delta_in_unit_interval_and_gamma_positive():
    g = heston_mc_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL,
                         n_steps=15, n_paths=3_000, seed=4)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0


def test_variance_vegas_positive():
    # More initial or long-run variance raises a call.
    g = heston_mc_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL,
                         n_steps=15, n_paths=3_000, seed=5)
    assert g["vega_v0"] > 0.0
    assert g["vega_theta"] > 0.0


def test_reproducible():
    kw = dict(n_steps=12, n_paths=1_500, seed=99)
    a = heston_mc_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, **kw)
    b = heston_mc_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, **kw)
    assert a["price"] == b["price"]
    assert a["delta"] == b["delta"]
