"""Tests for the Merton jump-diffusion implied-vol smile."""

import pytest

from quantforge import merton_smile


STRIKES = [80, 90, 100, 110, 120]


def test_zero_intensity_is_flat_at_sigma():
    sm = merton_smile(100, STRIKES, 1.0, 0.05, 0.2, lam=0.0, mu_j=0.0, sigma_j=0.0)
    for _, iv in sm:
        assert iv == pytest.approx(0.2, abs=1e-4)


def test_symmetric_jumps_give_a_smile():
    # Zero mean jump but positive jump vol fattens both tails -> wings above ATM.
    sm = merton_smile(100, STRIKES, 1.0, 0.05, 0.2, lam=1.0, mu_j=0.0, sigma_j=0.15)
    vols = [iv for _, iv in sm]
    atm = vols[len(vols) // 2]
    assert vols[0] > atm and vols[-1] > atm


def test_negative_mean_jump_gives_downward_skew():
    sm = merton_smile(100, STRIKES, 1.0, 0.05, 0.2, lam=1.0, mu_j=-0.2, sigma_j=0.1)
    vols = [iv for _, iv in sm]
    assert vols[0] > vols[-1]


def test_atm_vol_above_diffusion_vol():
    # Adding jumps raises the ATM implied vol above the pure diffusion vol.
    sm = merton_smile(100, [100], 1.0, 0.05, 0.2, lam=2.0, mu_j=0.0, sigma_j=0.2)
    assert sm[0][1] > 0.2


def test_all_vols_positive_and_sorted_axis():
    sm = merton_smile(100, STRIKES, 0.5, 0.03, 0.25, lam=1.5, mu_j=-0.1, sigma_j=0.15)
    ks = [k for k, _ in sm]
    assert ks == sorted(ks)
    assert all(iv > 0 for _, iv in sm)
