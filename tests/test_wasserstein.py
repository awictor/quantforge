"""Wasserstein-1 distance between two smile-implied densities."""

import math

import pytest

from quantforge import wasserstein_smiles as w1


S0, T, R = 100.0, 1.0, 0.05


def _flat(sig):
    return lambda K: sig


def test_identical_is_zero():
    assert w1(S0, T, R, _flat(0.2), _flat(0.2)) == pytest.approx(0.0, abs=1e-3)


def test_symmetric():
    a = w1(S0, T, R, _flat(0.2), _flat(0.3))
    b = w1(S0, T, R, _flat(0.3), _flat(0.2))
    assert a == pytest.approx(b, abs=1e-6)


def test_positive_for_different_vols():
    assert w1(S0, T, R, _flat(0.2), _flat(0.3)) > 0.0


def test_grows_with_vol_gap():
    near = w1(S0, T, R, _flat(0.2), _flat(0.22))
    far = w1(S0, T, R, _flat(0.2), _flat(0.35))
    assert far > near


def test_skew_vs_flat_positive():
    def down(K):
        return max(0.05, 0.2 + 0.15 * math.log(100.0 / K))
    assert w1(S0, T, R, down, _flat(0.2)) > 0.0


def test_triangle_inequality():
    # W1 is a metric: d(p, s) <= d(p, q) + d(q, s).
    p, q, s = _flat(0.18), _flat(0.24), _flat(0.30)
    d_ps = w1(S0, T, R, p, s)
    d_pq = w1(S0, T, R, p, q)
    d_qs = w1(S0, T, R, q, s)
    assert d_ps <= d_pq + d_qs + 1e-6
