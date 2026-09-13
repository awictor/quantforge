"""GPH long-memory estimator and fractional integration."""

import random

import pytest

from quantforge import gph_estimate, fractional_integrate
from quantforge.fracdiff import fractional_difference


def test_white_noise_near_zero_d():
    rng = random.Random(1)
    wn = [rng.gauss(0, 1) for _ in range(1024)]
    r = gph_estimate(wn)
    assert abs(r["d"]) < 3 * r["std_error"]     # not significantly different from 0


def test_std_error_positive_and_shrinks_with_m():
    rng = random.Random(2)
    x = [rng.gauss(0, 1) for _ in range(1024)]
    se_small = gph_estimate(x, m=20)["std_error"]
    se_large = gph_estimate(x, m=100)["std_error"]
    assert se_large < se_small           # more frequencies -> tighter


def test_integrate_then_difference_roundtrips():
    rng = random.Random(7)
    noise = [rng.gauss(0, 1) for _ in range(400)]
    back = fractional_difference(fractional_integrate(noise, 0.4), 0.4)
    assert max(abs(back[i] - noise[i]) for i in range(400)) < 1e-9


def test_integrate_d_zero_is_identity():
    rng = random.Random(3)
    noise = [rng.gauss(0, 1) for _ in range(50)]
    assert fractional_integrate(noise, 0.0) == noise


@pytest.mark.slow
def test_recovers_known_memory_parameter():
    for dtrue in (0.2, 0.4):
        ds = []
        for s in range(8):
            rng = random.Random(s)
            noise = [rng.gauss(0, 1) for _ in range(1024)]
            ds.append(gph_estimate(fractional_integrate(noise, dtrue))["d"])
        mean = sum(ds) / len(ds)
        assert abs(mean - dtrue) < 0.1     # GPH is slightly upward-biased


def test_validation():
    with pytest.raises(ValueError):
        gph_estimate([1.0, 2.0, 3.0])          # < 8 obs
    with pytest.raises(ValueError):
        gph_estimate([float(i) for i in range(64)], m=1)
    with pytest.raises(ValueError):
        fractional_integrate([], 0.3)
