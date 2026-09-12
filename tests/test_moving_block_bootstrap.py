"""Moving-block bootstrap CI for serially-correlated data."""

import random

import pytest

from quantforge import moving_block_bootstrap_ci, bootstrap_ci


def _ar1(n, phi, seed):
    rng = random.Random(seed)
    y = [0.0]
    for _ in range(n):
        y.append(phi * y[-1] + rng.gauss(0, 1))
    return y[1:]


def test_ci_brackets_sample_mean():
    ar = _ar1(2000, 0.7, 1)
    lo, mid, hi = moving_block_bootstrap_ci(ar, block=20)
    assert lo < mid < hi


def test_wider_than_iid_on_autocorrelated_mean():
    ar = _ar1(2000, 0.7, 1)
    lo, _, hi = moving_block_bootstrap_ci(ar, block=20)
    liid, _, hiid = bootstrap_ci(ar)
    assert (hi - lo) > (hiid - liid)          # iid under-covers correlated data


def test_block_one_matches_iid_width():
    ar = _ar1(2000, 0.7, 1)
    l1, _, h1 = moving_block_bootstrap_ci(ar, block=1)
    liid, _, hiid = bootstrap_ci(ar)
    assert abs((h1 - l1) - (hiid - liid)) / (hiid - liid) < 0.2


def test_iid_data_similar_to_iid_bootstrap():
    rng = random.Random(2)
    wn = [rng.gauss(0, 1) for _ in range(1000)]
    lw, _, hw = moving_block_bootstrap_ci(wn, block=20)
    liw, _, hiw = bootstrap_ci(wn)
    assert abs((hw - lw) - (hiw - liw)) / (hiw - liw) < 0.3


def test_validation():
    ar = _ar1(100, 0.5, 1)
    with pytest.raises(ValueError):
        moving_block_bootstrap_ci(ar, block=0)
    with pytest.raises(ValueError):
        moving_block_bootstrap_ci(ar, block=len(ar) + 1)
    with pytest.raises(ValueError):
        moving_block_bootstrap_ci([], block=1)
