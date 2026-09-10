"""Analytic single-slice / surface SVI Dupire local vol vs finite differences."""

import math

import pytest

from quantforge import (
    SVIParams,
    svi_local_variance,
    svi_surface_local_vol,
)
from quantforge.svi import _svi_derivs


def _slice_at(t):
    # Total variance scaling roughly with t (a simple linear term structure).
    return SVIParams(a=0.04 * t, b=0.4 * t, rho=-0.4, m=0.0, s=0.1)


SLICES = {0.5: _slice_at(0.5), 1.0: _slice_at(1.0), 2.0: _slice_at(2.0)}


def _w_kt(k, t):
    ts = sorted(SLICES)
    if t <= ts[0]:
        i = 0
    elif t >= ts[-1]:
        i = len(ts) - 2
    else:
        i = max(j for j in range(len(ts) - 1) if ts[j] <= t)
    t0, t1 = ts[i], ts[i + 1]
    lam = (t - t0) / (t1 - t0)
    return (1 - lam) * SLICES[t0].total_variance(k) + lam * SLICES[t1].total_variance(k)


def _dupire_fd(k, t):
    h, ht = 1e-4, 1e-4
    w = _w_kt(k, t)
    w_k = (_w_kt(k + h, t) - _w_kt(k - h, t)) / (2 * h)
    w_kk = (_w_kt(k + h, t) - 2 * w + _w_kt(k - h, t)) / h ** 2
    w_t = (_w_kt(k, t + ht) - _w_kt(k, t - ht)) / (2 * ht)
    denom = (1 - (k / w) * w_k
             + 0.25 * (-0.25 - 1 / w + k * k / (w * w)) * w_k * w_k
             + 0.5 * w_kk)
    return math.sqrt(w_t / denom)


@pytest.mark.parametrize("t", [0.75, 1.0, 1.5])
@pytest.mark.parametrize("k", [-0.15, 0.0, 0.15])
def test_surface_local_vol_matches_fd(t, k):
    assert svi_surface_local_vol(SLICES, k, t) == pytest.approx(_dupire_fd(k, t),
                                                                abs=1e-5)


def test_single_slice_uses_supplied_dw_dt():
    p = _slice_at(1.0)
    k = 0.1
    w, wp, wpp = _svi_derivs(p, k)
    dw_dt = 0.05
    denom = (1 - (k / w) * wp
             + 0.25 * (-0.25 - 1 / w + k * k / (w * w)) * wp * wp
             + 0.5 * wpp)
    assert svi_local_variance(p, k, dw_dt) == pytest.approx(dw_dt / denom)


def test_local_variance_positive_atm():
    assert svi_local_variance(_slice_at(1.0), 0.0, 0.04) > 0.0


def test_surface_rejects_out_of_range():
    with pytest.raises(ValueError):
        svi_surface_local_vol(SLICES, 0.0, 5.0)


def test_surface_needs_two_slices():
    with pytest.raises(ValueError):
        svi_surface_local_vol({1.0: _slice_at(1.0)}, 0.0, 1.0)


def test_negative_denominator_raises():
    # A wildly over-steep slice drives the Dupire denominator negative.
    bad = SVIParams(a=0.001, b=8.0, rho=-0.95, m=0.0, s=0.02)
    with pytest.raises(ValueError):
        svi_local_variance(bad, -0.3, 0.04)  # deep left wing: g < 0 there
