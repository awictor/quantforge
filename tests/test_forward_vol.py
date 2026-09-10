"""Tests for forward variance / forward vol on the term-structure surface."""

import math

import pytest

from quantforge import VolSurface, SurfaceSlice, SVIParams


def _flat(t, w):
    # A flat smile with constant total variance w.
    return SurfaceSlice(t=t, params=SVIParams(a=w, b=0.0, rho=0.0, m=0.0, s=0.1),
                        rmse=0.0)


def test_forward_variance_is_total_variance_increment():
    s = VolSurface([_flat(1.0, 0.04), _flat(2.0, 0.10)])
    # (0.10 - 0.04) / (2 - 1) = 0.06.
    assert s.forward_variance(0.0, 1.0, 2.0) == pytest.approx(0.06, abs=1e-12)


def test_forward_vol_is_sqrt_forward_variance():
    s = VolSurface([_flat(1.0, 0.04), _flat(2.0, 0.10)])
    assert s.forward_vol(0.0, 1.0, 2.0) == pytest.approx(math.sqrt(0.06), abs=1e-12)


def test_forward_vol_from_zero_equals_spot_vol():
    s = VolSurface([_flat(1.0, 0.04), _flat(2.0, 0.10)])
    assert s.forward_vol(0.0, 1e-9, 1.0) == pytest.approx(s.implied_vol(0.0, 1.0), abs=1e-4)


def test_forward_variance_additive():
    # Forward variances over adjacent windows reconstruct the total variance.
    s = VolSurface([_flat(1.0, 0.04), _flat(2.0, 0.10), _flat(3.0, 0.18)])
    fv1 = s.forward_variance(0.0, 1.0, 2.0)
    fv2 = s.forward_variance(0.0, 2.0, 3.0)
    total_23 = fv1 * 1.0 + fv2 * 1.0
    assert total_23 == pytest.approx(0.18 - 0.04, abs=1e-9)


def test_negative_forward_variance_raises():
    # Calendar-arbitrage surface: total variance falls with maturity.
    s = VolSurface([_flat(1.0, 0.10), _flat(2.0, 0.04)])
    with pytest.raises(ValueError):
        s.forward_vol(0.0, 1.0, 2.0)


def test_requires_ordered_times():
    s = VolSurface([_flat(1.0, 0.04), _flat(2.0, 0.10)])
    with pytest.raises(ValueError):
        s.forward_variance(0.0, 2.0, 1.0)
