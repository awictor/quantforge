"""Constant Proportion Portfolio Insurance."""

import math

import pytest

from quantforge import discounted_floor, cushion, risky_exposure, cppi_path


def test_discounted_floor():
    assert discounted_floor(90, 0.03, 2) == pytest.approx(90 * math.exp(-0.06))


def test_cushion():
    assert cushion(100, 85) == 15
    assert cushion(80, 85) == 0


def test_risky_exposure_multiplier_and_caps():
    assert risky_exposure(100, 85, 3) == 45      # 3 * 15
    assert risky_exposure(100, 50, 5) == 100     # capped at wealth
    assert risky_exposure(84, 85, 3) == 0        # cushion zero


def test_higher_multiplier_more_exposure():
    assert risky_exposure(100, 85, 4) > risky_exposure(100, 85, 2)


def test_path_grows_with_positive_returns():
    p = cppi_path(100, 90, 3, [0.05] * 4, 0.02, 0.25)
    assert p[-1] > 100


def test_path_protected_in_crash():
    crash = cppi_path(100, 90, 3, [-0.30] * 4, 0.02, 0.25)
    assert crash[-1] >= 85   # near the guaranteed floor


def test_validation():
    with pytest.raises(ValueError):
        risky_exposure(-1, 85, 3)
    with pytest.raises(ValueError):
        cppi_path(-1, 90, 3, [0.01], 0.02, 0.25)
