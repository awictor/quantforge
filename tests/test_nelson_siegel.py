"""Nelson-Siegel and Svensson parametric yield curves."""

import math

import pytest

from quantforge import (
    nelson_siegel_zero, svensson_zero, nelson_siegel_discount,
    nelson_siegel_forward,
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


def test_validation():
    with pytest.raises(ValueError):
        nelson_siegel_zero(5, B0, B1, B2, 0)
    with pytest.raises(ValueError):
        svensson_zero(5, B0, B1, B2, 0.0, TAU, 0)
