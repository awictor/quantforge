"""Tests for the Heston implied-vol smile extraction."""

import math

import pytest

from quantforge import heston_smile


STRIKES = [80, 90, 100, 110, 120]


def test_zero_volvol_is_flat_at_sqrt_v0():
    sm = heston_smile(100, STRIKES, 1.0, 0.0, v0=0.04, kappa=2.0, theta=0.04,
                      xi=1e-5, rho=0.0)
    for _, iv in sm:
        assert iv == pytest.approx(0.2, abs=1e-3)


def test_negative_rho_gives_downward_skew():
    sm = heston_smile(100, STRIKES, 1.0, 0.0, v0=0.04, kappa=2.0, theta=0.04,
                      xi=0.5, rho=-0.7)
    vols = [iv for _, iv in sm]
    # Monotone decreasing in strike (downward equity skew).
    assert all(vols[i] > vols[i + 1] for i in range(len(vols) - 1))


def test_positive_rho_gives_upward_skew():
    sm = heston_smile(100, STRIKES, 1.0, 0.0, v0=0.04, kappa=2.0, theta=0.04,
                      xi=0.5, rho=0.7)
    vols = [iv for _, iv in sm]
    assert all(vols[i] < vols[i + 1] for i in range(len(vols) - 1))


def test_log_moneyness_axis_sorted_and_centered():
    sm = heston_smile(100, STRIKES, 1.0, 0.0, v0=0.04, kappa=2.0, theta=0.04,
                      xi=0.4, rho=-0.3)
    ks = [k for k, _ in sm]
    assert ks == sorted(ks)
    # ATM strike 100 with r=0 -> forward 100 -> k=0.
    atm = [iv for k, iv in sm if abs(k) < 1e-9]
    assert atm  # ATM point present


def test_all_vols_positive():
    sm = heston_smile(100, STRIKES, 0.5, 0.03, v0=0.05, kappa=1.5, theta=0.06,
                      xi=0.6, rho=-0.5)
    assert all(iv > 0 for _, iv in sm)
