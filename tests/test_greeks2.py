"""Second-order Greeks verified against finite differences.

Each analytic higher-order Greek is checked against a central finite difference
of the corresponding first-order Greek (or price), for both calls and puts and
across several strikes. This catches any algebra error in the closed forms.
"""

import pytest

from quantforge import (
    delta, vega, gamma, OptionType,
    vanna, vomma, volga, charm, veta, speed, zomma, color,
)

CASES = [
    dict(S=100, K=90, t=0.75, r=0.03, sigma=0.25),
    dict(S=100, K=100, t=0.75, r=0.03, sigma=0.25),
    dict(S=100, K=115, t=0.40, r=0.05, sigma=0.35),
]


@pytest.mark.parametrize("p", CASES)
def test_vanna_matches_fd_of_vega_in_spot(p):
    h = 1e-3
    fd = (vega(p["S"] + h, p["K"], p["t"], p["r"], p["sigma"])
          - vega(p["S"] - h, p["K"], p["t"], p["r"], p["sigma"])) / (2 * h)
    assert vanna(**p) == pytest.approx(fd, abs=1e-4)


@pytest.mark.parametrize("p", CASES)
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_vanna_matches_fd_of_delta_in_vol(p, ot):
    h = 1e-5
    fd = (delta(p["S"], p["K"], p["t"], p["r"], p["sigma"] + h, ot)
          - delta(p["S"], p["K"], p["t"], p["r"], p["sigma"] - h, ot)) / (2 * h)
    assert vanna(**p) == pytest.approx(fd, abs=1e-3)


@pytest.mark.parametrize("p", CASES)
def test_vomma_matches_fd_of_vega_in_vol(p):
    h = 1e-5
    fd = (vega(p["S"], p["K"], p["t"], p["r"], p["sigma"] + h)
          - vega(p["S"], p["K"], p["t"], p["r"], p["sigma"] - h)) / (2 * h)
    assert vomma(**p) == pytest.approx(fd, abs=1e-3)
    assert volga(**p) == vomma(**p)


@pytest.mark.parametrize("p", CASES)
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_charm_matches_fd_of_delta_in_time(p, ot):
    h = 1e-5
    # Calendar charm = -d(delta)/d(t_expiry).
    fd = -(delta(p["S"], p["K"], p["t"] + h, p["r"], p["sigma"], ot)
           - delta(p["S"], p["K"], p["t"] - h, p["r"], p["sigma"], ot)) / (2 * h)
    assert charm(**p, option_type=ot) == pytest.approx(fd, abs=1e-3)


@pytest.mark.parametrize("p", CASES)
def test_veta_matches_fd_of_vega_in_time(p):
    h = 1e-5
    fd = -(vega(p["S"], p["K"], p["t"] + h, p["r"], p["sigma"])
           - vega(p["S"], p["K"], p["t"] - h, p["r"], p["sigma"])) / (2 * h)
    assert veta(**p) == pytest.approx(fd, abs=1e-2)


@pytest.mark.parametrize("p", CASES)
def test_speed_matches_fd_of_gamma_in_spot(p):
    h = 1e-2
    fd = (gamma(p["S"] + h, p["K"], p["t"], p["r"], p["sigma"])
          - gamma(p["S"] - h, p["K"], p["t"], p["r"], p["sigma"])) / (2 * h)
    assert speed(**p) == pytest.approx(fd, abs=1e-5)


@pytest.mark.parametrize("p", CASES)
def test_zomma_matches_fd_of_gamma_in_vol(p):
    h = 1e-5
    fd = (gamma(p["S"], p["K"], p["t"], p["r"], p["sigma"] + h)
          - gamma(p["S"], p["K"], p["t"], p["r"], p["sigma"] - h)) / (2 * h)
    assert zomma(**p) == pytest.approx(fd, abs=1e-4)


@pytest.mark.parametrize("p", CASES)
def test_color_matches_fd_of_gamma_in_time(p):
    h = 1e-5
    fd = -(gamma(p["S"], p["K"], p["t"] + h, p["r"], p["sigma"])
           - gamma(p["S"], p["K"], p["t"] - h, p["r"], p["sigma"])) / (2 * h)
    assert color(**p) == pytest.approx(fd, abs=1e-3)
