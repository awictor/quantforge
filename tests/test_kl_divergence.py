"""Kullback-Leibler divergence between two smile-implied densities."""

import math

import pytest

from quantforge import kl_divergence_smiles as kl


S0, T, R = 100.0, 1.0, 0.05


def _flat(sig):
    return lambda K: sig


def test_identical_smiles_zero():
    assert kl(S0, T, R, _flat(0.2), _flat(0.2)) == pytest.approx(0.0, abs=1e-4)


def test_different_vols_positive():
    assert kl(S0, T, R, _flat(0.2), _flat(0.3)) > 0.0


def test_asymmetric():
    a = kl(S0, T, R, _flat(0.2), _flat(0.3))
    b = kl(S0, T, R, _flat(0.3), _flat(0.2))
    assert abs(a - b) > 1e-3


def test_non_negative():
    for p, q in [(0.15, 0.25), (0.3, 0.2), (0.2, 0.22)]:
        assert kl(S0, T, R, _flat(p), _flat(q)) >= -1e-6


def test_skew_vs_flat_positive():
    def down(K):
        return max(0.05, 0.2 + 0.15 * math.log(100.0 / K))
    assert kl(S0, T, R, down, _flat(0.2)) > 0.0


def test_larger_vol_gap_larger_divergence():
    near = kl(S0, T, R, _flat(0.2), _flat(0.22))
    far = kl(S0, T, R, _flat(0.2), _flat(0.35))
    assert far > near
