"""Analytic SSVI Dupire local variance vs a finite-difference cross-check."""

import math

import pytest

from quantforge import (
    SSVIParams,
    ssvi_total_variance,
    ssvi_local_variance,
    ssvi_local_vol_from_params,
)


RHO, ETA, GAMMA = -0.4, 1.0, 0.5
C = 0.04  # theta(t) = C * t, so theta'(t) = C


def _w(k, t):
    return ssvi_total_variance(k, C * t, RHO, ETA, GAMMA)


def _dupire_fd(k, t):
    h, ht = 1e-4, 1e-4
    w = _w(k, t)
    w_k = (_w(k + h, t) - _w(k - h, t)) / (2 * h)
    w_kk = (_w(k + h, t) - 2 * w + _w(k - h, t)) / h ** 2
    w_t = (_w(k, t + ht) - _w(k, t - ht)) / (2 * ht)
    denom = (1 - (k / w) * w_k
             + 0.25 * (-0.25 - 1 / w + k * k / (w * w)) * w_k * w_k
             + 0.5 * w_kk)
    return w_t / denom


@pytest.mark.parametrize("t", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("k", [-0.2, 0.0, 0.2])
def test_analytic_local_variance_matches_fd(t, k):
    an = ssvi_local_variance(k, t, C * t, C, RHO, ETA, GAMMA)
    assert an == pytest.approx(_dupire_fd(k, t), abs=1e-6)


def test_atm_local_variance_positive():
    lv = ssvi_local_variance(0.0, 1.0, C * 1.0, C, RHO, ETA, GAMMA)
    assert lv > 0.0


def test_from_params_interpolates_theta():
    params = SSVIParams(rho=RHO, eta=ETA, gamma=GAMMA,
                        thetas={0.5: C * 0.5, 1.0: C * 1.0, 2.0: C * 2.0})
    # At a pillar it should match the direct analytic call.
    direct = math.sqrt(ssvi_local_variance(0.1, 1.0, C * 1.0, C, RHO, ETA, GAMMA))
    assert ssvi_local_vol_from_params(params, 0.1, 1.0) == pytest.approx(direct,
                                                                         abs=1e-9)


def test_from_params_rejects_out_of_range():
    params = SSVIParams(rho=RHO, eta=ETA, gamma=GAMMA,
                        thetas={0.5: 0.02, 1.0: 0.04})
    with pytest.raises(ValueError):
        ssvi_local_vol_from_params(params, 0.0, 3.0)


def test_negative_denominator_raises():
    # A pathological slice with a huge skew drives the Dupire denominator
    # negative (butterfly arbitrage); the analytic form should flag it.
    with pytest.raises(ValueError):
        ssvi_local_variance(0.5, 0.5, 0.5, 0.04, -0.9, 40.0, 0.2)
