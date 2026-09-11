"""Greeks of range binary (analytic delta/gamma) and supershare (FD)."""

import pytest

from quantforge import (
    range_binary, range_binary_greeks, supershare, supershare_greeks,
)


CASES = [
    (100.0, 95.0, 110.0, 1.0, 0.05, 0.2),
    (100.0, 90.0, 105.0, 0.5, 0.03, 0.3),
    (80.0, 70.0, 120.0, 2.0, 0.02, 0.25),
]


@pytest.mark.parametrize("S,Kl,Kh,t,r,sig", CASES)
def test_range_binary_delta_matches_fd(S, Kl, Kh, t, r, sig):
    g = range_binary_greeks(S, Kl, Kh, t, r, sig)
    h = 1e-3 * S
    fd = (range_binary(S + h, Kl, Kh, t, r, sig)
          - range_binary(S - h, Kl, Kh, t, r, sig)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-6)


@pytest.mark.parametrize("S,Kl,Kh,t,r,sig", CASES)
def test_range_binary_gamma_matches_fd(S, Kl, Kh, t, r, sig):
    g = range_binary_greeks(S, Kl, Kh, t, r, sig)
    h = 1e-3 * S
    fd = (range_binary(S + h, Kl, Kh, t, r, sig)
          - 2 * range_binary(S, Kl, Kh, t, r, sig)
          + range_binary(S - h, Kl, Kh, t, r, sig)) / (h * h)
    assert g["gamma"] == pytest.approx(fd, abs=1e-5)


def test_range_binary_price_field_matches_pricer():
    g = range_binary_greeks(100.0, 95.0, 110.0, 1.0, 0.05, 0.2)
    assert g["price"] == pytest.approx(
        range_binary(100.0, 95.0, 110.0, 1.0, 0.05, 0.2), abs=1e-12)


def test_range_binary_theta_positive_when_centered():
    # A corridor centered on the forward gains value as expiry nears: the
    # terminal mass concentrates inside, so P(in corridor) -> 1 and the price
    # rises toward the discounted cash. Calendar theta = -dV/dt > 0.
    fwd = 100.0 * 2.718281828 ** (0.05 * 1.0)
    g = range_binary_greeks(100.0, fwd - 8.0, fwd + 8.0, 1.0, 0.05, 0.2)
    assert g["theta"] > 0.0


def test_supershare_price_field_and_delta_matches_fd():
    S, Kl, Kh, t, r, sig = 100.0, 95.0, 110.0, 1.0, 0.05, 0.2
    g = supershare_greeks(S, Kl, Kh, t, r, sig)
    assert g["price"] == pytest.approx(supershare(S, Kl, Kh, t, r, sig), abs=1e-12)
    h = 1e-4 * S
    fd = (supershare(S + h, Kl, Kh, t, r, sig)
          - supershare(S - h, Kl, Kh, t, r, sig)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-5)


def test_validation():
    with pytest.raises(ValueError):
        range_binary_greeks(100.0, 110.0, 95.0, 1.0, 0.05, 0.2)
    with pytest.raises(ValueError):
        supershare_greeks(100.0, 110.0, 95.0, 1.0, 0.05, 0.2)
