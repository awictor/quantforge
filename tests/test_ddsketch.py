"""Tests for DDSketch (relative-error quantile sketch).

Estimates are cross-checked against the exact sorted order statistics, verifying the
relative-error guarantee rather than exact equality.
"""

import math
import random

import pytest

from quantforge.ddsketch import DDSketch


def _lognormal(seed=3, n=50000, sigma=2.0):
    rng = random.Random(seed)
    return [math.exp(rng.gauss(0, sigma)) for _ in range(n)]


def _true_quantile(sorted_xs, q):
    rank = q * (len(sorted_xs) - 1)
    return sorted_xs[int(math.floor(rank))]


def test_relative_error_within_alpha():
    alpha = 0.01
    xs = _lognormal()
    s = DDSketch(alpha)
    for x in xs:
        s.add(x)
    sx = sorted(xs)
    for q in (0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999):
        est = s.quantile(q)
        true = _true_quantile(sx, q)
        assert abs(est - true) / true <= alpha + 1e-9


def test_relative_error_larger_alpha():
    alpha = 0.05
    xs = _lognormal(seed=9, n=20000)
    s = DDSketch(alpha)
    for x in xs:
        s.add(x)
    sx = sorted(xs)
    for q in (0.05, 0.5, 0.95):
        est = s.quantile(q)
        true = _true_quantile(sx, q)
        assert abs(est - true) / true <= alpha + 1e-9


def test_min_and_max_are_exact():
    xs = _lognormal(n=10000)
    s = DDSketch(0.01)
    for x in xs:
        s.add(x)
    sx = sorted(xs)
    assert s.min() == sx[0]   # min()/max() are exact accessors
    assert s.max() == sx[-1]
    # extreme quantiles are relative-error estimates clamped to the observed range
    assert sx[0] <= s.quantile(0.0)
    assert abs(s.quantile(0.0) - sx[0]) / sx[0] <= 0.01
    assert s.quantile(1.0) <= sx[-1]
    assert abs(s.quantile(1.0) - sx[-1]) / sx[-1] <= 0.01


def test_merge_matches_whole():
    alpha = 0.01
    xs = _lognormal(seed=5, n=40000)
    whole = DDSketch(alpha)
    for x in xs:
        whole.add(x)
    a = DDSketch(alpha)
    b = DDSketch(alpha)
    for i, x in enumerate(xs):
        (a if i % 2 else b).add(x)
    m = a + b
    assert m.n == whole.n
    for q in (0.1, 0.5, 0.9, 0.99):
        assert m.quantile(q) == pytest.approx(whole.quantile(q), abs=1e-9)


def test_merge_in_place():
    a = DDSketch(0.02)
    b = DDSketch(0.02)
    for x in (1.0, 2.0, 3.0):
        a.add(x)
    for x in (4.0, 5.0, 6.0):
        b.add(x)
    a.merge(b)
    assert a.n == 6
    assert a.min() == 1.0
    assert a.max() == 6.0


def test_memory_sublinear_in_n():
    xs = _lognormal(n=50000)
    s = DDSketch(0.01)
    for x in xs:
        s.add(x)
    # buckets grow with the log-range, far below the sample size
    assert s.num_buckets() < len(xs) // 10


def test_weighted_add_matches_repeated():
    a = DDSketch(0.01)
    a.add(10.0, weight=5)
    b = DDSketch(0.01)
    for _ in range(5):
        b.add(10.0)
    assert a.n == b.n == 5
    assert a.quantile(0.5) == pytest.approx(b.quantile(0.5))


def test_single_value_quantiles():
    s = DDSketch(0.01)
    s.add(42.0)
    assert abs(s.quantile(0.5) - 42.0) / 42.0 <= 0.01
    assert s.quantile(0.0) == 42.0
    assert s.quantile(1.0) == 42.0


def test_uniform_data_median():
    s = DDSketch(0.01)
    xs = [float(i) for i in range(1, 10001)]
    for x in xs:
        s.add(x)
    est = s.quantile(0.5)
    true = xs[int(0.5 * (len(xs) - 1))]
    assert abs(est - true) / true <= 0.01


def test_alpha_out_of_range_raises():
    with pytest.raises(ValueError):
        DDSketch(0.0)
    with pytest.raises(ValueError):
        DDSketch(1.0)


def test_nonpositive_value_raises():
    s = DDSketch(0.01)
    with pytest.raises(ValueError):
        s.add(0.0)
    with pytest.raises(ValueError):
        s.add(-1.0)


def test_nonpositive_weight_raises():
    with pytest.raises(ValueError):
        DDSketch(0.01).add(1.0, weight=0)


def test_quantile_out_of_range_raises():
    s = DDSketch(0.01)
    s.add(1.0)
    with pytest.raises(ValueError):
        s.quantile(-0.1)
    with pytest.raises(ValueError):
        s.quantile(1.1)


def test_empty_sketch_raises():
    s = DDSketch(0.01)
    with pytest.raises(ValueError):
        s.quantile(0.5)
    with pytest.raises(ValueError):
        s.min()
    with pytest.raises(ValueError):
        s.max()


def test_merge_alpha_mismatch_raises():
    with pytest.raises(ValueError):
        DDSketch(0.01) + DDSketch(0.02)


def test_merge_wrong_type_raises():
    with pytest.raises(TypeError):
        DDSketch(0.01).merge([1, 2, 3])


def test_add_returns_self():
    s = DDSketch(0.01)
    assert s.add(1.0) is s
