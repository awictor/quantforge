"""Robust scale estimators: Qn, Sn, biweight midvariance."""

import random
import statistics

import pytest

from quantforge import qn_scale, sn_scale, biweight_midvariance


@pytest.mark.slow
def test_consistent_with_sigma_on_gaussian():
    rng = random.Random(1)
    data = [rng.gauss(0, 2.0) for _ in range(3000)]
    for est in (qn_scale, sn_scale, biweight_midvariance):
        assert abs(est(data) - 2.0) < 0.2


def test_robust_to_contamination():
    rng = random.Random(2)
    clean = [rng.gauss(0, 1.0) for _ in range(800)]
    cont = clean + [1000.0] * 200                       # 20% gross outliers
    # standard deviation is destroyed; robust estimators barely move
    assert statistics.pstdev(cont) > 100
    for est in (qn_scale, sn_scale):
        assert est(cont) < 2.0                          # still near the clean scale


def test_sn_known_value():
    # Sn of 1..10 (no finite-sample correction) matches robustbase's asymptotic form
    assert abs(sn_scale(list(range(1, 11))) - 3.5778) < 1e-3


@pytest.mark.slow
def test_scale_tracks_sigma():
    rng = random.Random(3)
    d = [rng.gauss(0, 5.0) for _ in range(5000)]
    assert abs(qn_scale(d) - 5.0) < 0.4
    assert abs(sn_scale(d) - 5.0) < 0.4


def test_validation():
    for est in (qn_scale, sn_scale, biweight_midvariance):
        with pytest.raises(ValueError):
            est([1.0])
