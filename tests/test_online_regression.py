"""Tests for RunningRegression (online simple OLS)."""

import random
import statistics

import pytest

from quantforge.online_regression import RunningRegression


def test_matches_statistics_linear_regression():
    random.seed(7)
    xs = [random.uniform(-5, 5) for _ in range(500)]
    ys = [3.2 * x - 1.7 + random.gauss(0, 0.5) for x in xs]
    r = RunningRegression(xs, ys)
    lr = statistics.linear_regression(xs, ys)
    assert r.slope() == pytest.approx(lr.slope, abs=1e-9)
    assert r.intercept() == pytest.approx(lr.intercept, abs=1e-9)


def test_r_squared_equals_correlation_squared():
    random.seed(11)
    xs = [random.uniform(0, 10) for _ in range(300)]
    ys = [1.5 * x + 4 + random.gauss(0, 2) for x in xs]
    r = RunningRegression(xs, ys)
    cc = statistics.correlation(xs, ys)
    assert r.correlation() == pytest.approx(cc, abs=1e-9)
    assert r.r_squared() == pytest.approx(cc * cc, abs=1e-9)


def test_perfect_line_recovered_exactly():
    r = RunningRegression([1, 2, 3, 4], [2 * x + 5 for x in (1, 2, 3, 4)])
    assert r.slope() == pytest.approx(2.0, abs=1e-12)
    assert r.intercept() == pytest.approx(5.0, abs=1e-12)
    assert r.r_squared() == pytest.approx(1.0, abs=1e-12)
    assert r.predict(10) == pytest.approx(25.0, abs=1e-9)


def test_incremental_update_matches_batch():
    random.seed(3)
    xs = [random.uniform(-2, 2) for _ in range(120)]
    ys = [-0.8 * x + 2.0 + random.gauss(0, 0.3) for x in xs]
    batch = RunningRegression(xs, ys)
    inc = RunningRegression()
    for x, y in zip(xs, ys):
        inc.update(x, y)
    assert inc.slope() == pytest.approx(batch.slope(), abs=1e-12)
    assert inc.intercept() == pytest.approx(batch.intercept(), abs=1e-12)


def test_merge_matches_whole_sample():
    random.seed(5)
    xs = [random.uniform(-5, 5) for _ in range(400)]
    ys = [2.5 * x - 0.5 + random.gauss(0, 1) for x in xs]
    whole = RunningRegression(xs, ys)
    a = RunningRegression(xs[:137], ys[:137])
    b = RunningRegression(xs[137:], ys[137:])
    m = a + b
    assert m.slope() == pytest.approx(whole.slope(), abs=1e-9)
    assert m.intercept() == pytest.approx(whole.intercept(), abs=1e-9)
    assert m.r_squared() == pytest.approx(whole.r_squared(), abs=1e-9)


def test_merge_with_empty_returns_other():
    r = RunningRegression([1, 2, 3], [2, 4, 6])
    e = RunningRegression()
    assert (e + r).slope() == pytest.approx(r.slope())
    assert (r + e).slope() == pytest.approx(r.slope())
    # merge result is an independent copy
    assert (e + r) is not r


def test_negative_slope():
    r = RunningRegression([0, 1, 2, 3], [10, 8, 6, 4])
    assert r.slope() == pytest.approx(-2.0)
    assert r.intercept() == pytest.approx(10.0)
    assert r.correlation() == pytest.approx(-1.0)


def test_slope_needs_two_points():
    with pytest.raises(ValueError):
        RunningRegression([1], [1]).slope()


def test_correlation_needs_two_points():
    with pytest.raises(ValueError):
        RunningRegression([1], [1]).correlation()


def test_zero_variance_x_raises():
    with pytest.raises(ValueError):
        RunningRegression([2, 2, 2], [1, 2, 3]).slope()


def test_zero_variance_correlation_raises():
    with pytest.raises(ValueError):
        RunningRegression([2, 2, 2], [1, 2, 3]).correlation()


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        RunningRegression([1, 2], [1])


def test_update_returns_self_for_chaining():
    r = RunningRegression()
    assert r.update(1, 2) is r


def test_correlation_clamped_to_unit_interval():
    r = RunningRegression([1, 2, 3, 4, 5], [2, 4, 6, 8, 10])
    assert -1.0 <= r.correlation() <= 1.0
    assert r.correlation() == pytest.approx(1.0)
