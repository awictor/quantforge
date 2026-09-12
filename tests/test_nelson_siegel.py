"""Nelson-Siegel and Svensson parametric yield curves."""

import math

import pytest

from quantforge import (
    nelson_siegel_zero, svensson_zero, nelson_siegel_discount,
    nelson_siegel_forward, fit_nelson_siegel,
)


B0, B1, B2, TAU = 0.04, -0.02, 0.01, 2.0


def test_short_rate_limit():
    assert nelson_siegel_zero(0, B0, B1, B2, TAU) == pytest.approx(B0 + B1)
    assert nelson_siegel_zero(1e-8, B0, B1, B2, TAU) == pytest.approx(B0 + B1, abs=1e-4)


def test_long_level_limit():
    assert nelson_siegel_zero(500, B0, B1, B2, TAU) == pytest.approx(B0, abs=1e-3)


def test_forward_limits():
    assert nelson_siegel_forward(0, B0, B1, B2, TAU) == pytest.approx(B0 + B1)
    assert nelson_siegel_forward(500, B0, B1, B2, TAU) == pytest.approx(B0, abs=1e-6)


def test_discount_from_zero():
    assert nelson_siegel_discount(5, B0, B1, B2, TAU) == pytest.approx(
        math.exp(-nelson_siegel_zero(5, B0, B1, B2, TAU) * 5))


def test_discount_decreasing():
    dfs = [nelson_siegel_discount(t, B0, B1, B2, TAU) for t in (1, 2, 5, 10, 20)]
    assert all(dfs[i] > dfs[i + 1] for i in range(len(dfs) - 1))


def test_svensson_reduces_to_ns():
    for t in (0.5, 2, 5, 10):
        assert svensson_zero(t, B0, B1, B2, 0.0, TAU, 3.0) == pytest.approx(
            nelson_siegel_zero(t, B0, B1, B2, TAU), abs=1e-12)


def test_svensson_second_hump_changes_curve():
    assert abs(svensson_zero(3, B0, B1, B2, 0.02, TAU, 5.0)
               - nelson_siegel_zero(3, B0, B1, B2, TAU)) > 1e-4


MATS = [0.5, 1, 2, 3, 5, 7, 10, 20, 30]


def test_fit_recovers_noiseless_parameters():
    zs = [nelson_siegel_zero(t, B0, B1, B2, TAU) for t in MATS]
    b0, b1, b2, tau = fit_nelson_siegel(MATS, zs)
    assert b0 == pytest.approx(B0, abs=1e-6)
    assert b1 == pytest.approx(B1, abs=1e-6)
    assert b2 == pytest.approx(B2, abs=1e-6)
    assert tau == pytest.approx(TAU, abs=1e-6)


def test_fit_reprices_curve():
    zs = [nelson_siegel_zero(t, B0, B1, B2, TAU) for t in MATS]
    b0, b1, b2, tau = fit_nelson_siegel(MATS, zs)
    assert all(nelson_siegel_zero(t, b0, b1, b2, tau) == pytest.approx(zs[i], abs=1e-6)
               for i, t in enumerate(MATS))


def test_fit_noisy_small_residual():
    zs = [nelson_siegel_zero(t, B0, B1, B2, TAU) for t in MATS]
    zn = [zs[i] + (0.0002 if i % 2 else -0.0002) for i in range(len(zs))]
    b0, b1, b2, tau = fit_nelson_siegel(MATS, zn)
    sse = sum((nelson_siegel_zero(t, b0, b1, b2, tau) - zn[i]) ** 2
              for i, t in enumerate(MATS))
    assert sse < 1e-5


def test_fit_validation():
    with pytest.raises(ValueError):
        fit_nelson_siegel([1, 2], [0.02, 0.03])


def test_validation():
    with pytest.raises(ValueError):
        nelson_siegel_zero(5, B0, B1, B2, 0)
    with pytest.raises(ValueError):
        svensson_zero(5, B0, B1, B2, 0.0, TAU, 0)
