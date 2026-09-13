"""MinRV, MedRV jump-robust variance and realized quarticity."""

import math
import random

import pytest

from quantforge import (min_realized_variance, med_realized_variance,
                        realized_quarticity, realized_variance_from_returns,
                        bipower_variation)


def _path(sigma, n, dt, jump=0.0, seed=0):
    rng = random.Random(seed)
    r = [sigma * math.sqrt(dt) * rng.gauss(0, 1) for _ in range(n)]
    if jump:
        r[n // 2] += jump
    return r


def test_no_jump_agrees_with_rv():
    sigma, n = 0.01, 2000
    iv = sigma * sigma
    mn = [min_realized_variance(_path(sigma, n, 1.0 / n, 0.0, s)) for s in range(40)]
    md = [med_realized_variance(_path(sigma, n, 1.0 / n, 0.0, s)) for s in range(40)]
    assert abs(sum(mn) / 40 - iv) < 0.1 * iv
    assert abs(sum(md) / 40 - iv) < 0.1 * iv


def test_jump_robustness():
    sigma, n, jump = 0.01, 2000, 0.05
    rv, mn, md = [], [], []
    for s in range(40):
        r = _path(sigma, n, 1.0 / n, jump, s)
        rv.append(realized_variance_from_returns(r))
        mn.append(min_realized_variance(r))
        md.append(med_realized_variance(r))
    mean_rv = sum(rv) / 40
    assert sum(mn) / 40 < 0.5 * mean_rv       # jump discarded
    assert sum(md) / 40 < 0.5 * mean_rv


def test_minrv_medrv_close_to_bipower_no_jump():
    r = _path(0.01, 3000, 1.0 / 3000, 0.0, 3)
    bv = bipower_variation(r)
    assert abs(min_realized_variance(r) - bv) < 0.15 * bv
    assert abs(med_realized_variance(r) - bv) < 0.15 * bv


def test_quarticity_positive_and_consistent():
    sigma, n = 0.01, 3000
    rq = [realized_quarticity(_path(sigma, n, 1.0 / n, 0.0, s)) for s in range(40)]
    mean_rq = sum(rq) / 40
    assert mean_rq > 0.0
    # Integrated quarticity of the path is sigma^4.
    assert abs(mean_rq - sigma ** 4) < 0.2 * sigma ** 4


def test_medrv_uses_median_of_three():
    # A single huge middle return is dropped by the median.
    base = [0.01, 0.01, 0.01, 0.01, 0.01]
    spiked = list(base)
    spiked[2] = 10.0
    assert abs(med_realized_variance(base) - med_realized_variance(spiked)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        min_realized_variance([0.01])
    with pytest.raises(ValueError):
        med_realized_variance([0.01, 0.02])
    with pytest.raises(ValueError):
        realized_quarticity([])
