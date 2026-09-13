"""Brier decomposition, reliability curve, expected calibration error."""

import random

import pytest

from quantforge import (
    brier_decomposition,
    reliability_curve,
    expected_calibration_error,
    brier_score,
)


def test_decomposition_identity_exact():
    # With grouping by identical value, reliability - resolution + uncertainty == Brier.
    rng = random.Random(4)
    for _ in range(300):
        n = rng.randint(2, 40)
        levels = [round(rng.random(), 2) for _ in range(rng.randint(2, 6))]
        f = [rng.choice(levels) for _ in range(n)]
        o = [1 if rng.random() < fi else 0 for fi in f]
        d = brier_decomposition(f, o)
        assert abs(d["brier"] - brier_score(o, f)) < 1e-9


def test_perfect_forecast_zero_reliability_and_brier():
    d = brier_decomposition([0.0, 0.0, 1.0, 1.0], [0, 0, 1, 1])
    assert abs(d["reliability"]) < 1e-12
    assert abs(d["brier"]) < 1e-12


def test_constant_at_base_rate_zero_resolution():
    o = [1, 0, 1, 0, 1, 0]
    d = brier_decomposition([0.5] * 6, o)
    assert abs(d["resolution"]) < 1e-12
    assert abs(d["reliability"]) < 1e-12
    assert abs(d["uncertainty"] - 0.25) < 1e-12
    assert abs(d["base_rate"] - 0.5) < 1e-12


def test_reliability_curve_on_diagonal_when_calibrated():
    rng = random.Random(1)
    f = [rng.random() for _ in range(20000)]
    o = [1 if rng.random() < fi else 0 for fi in f]
    curve = reliability_curve(f, o, n_bins=5)
    assert len(curve) == 5
    for mean_f, obs, cnt in curve:
        assert abs(mean_f - obs) < 0.03


def test_ece_calibrated_vs_miscalibrated():
    rng = random.Random(1)
    f = [rng.random() for _ in range(20000)]
    o = [1 if rng.random() < fi else 0 for fi in f]
    assert expected_calibration_error(f, o, n_bins=5) < 0.02

    f2 = [0.9] * 1000
    o2 = [1 if rng.random() < 0.5 else 0 for _ in range(1000)]
    assert expected_calibration_error(f2, o2, n_bins=5) > 0.3


def test_uncertainty_is_variance_of_outcome():
    o = [1, 1, 1, 0]  # base 0.75 -> 0.75*0.25 = 0.1875
    d = brier_decomposition([0.5] * 4, o)
    assert abs(d["uncertainty"] - 0.1875) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        brier_decomposition([0.5, 0.5], [1])          # length mismatch
    with pytest.raises(ValueError):
        brier_decomposition([1.5], [1])               # forecast out of range
    with pytest.raises(ValueError):
        brier_decomposition([0.5], [2])               # outcome not 0/1
