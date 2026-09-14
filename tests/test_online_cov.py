"""Online (single-pass) covariance and correlation."""

import random
import statistics

import pytest

from quantforge import RunningCovariance


def _batch_cov(xs, ys, ddof=1):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    return sum((xs[i] - mx) * (ys[i] - my) for i in range(n)) / (n - ddof)


def _batch_corr(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / (sxx * syy) ** 0.5


def test_vs_batch():
    rng = random.Random(1)
    for _ in range(500):
        n = rng.randint(2, 50)
        xs = [rng.uniform(-10, 10) for _ in range(n)]
        ys = [rng.uniform(-10, 10) for _ in range(n)]
        rc = RunningCovariance(xs, ys)
        assert abs(rc.covariance() - _batch_cov(xs, ys)) < 1e-7
        assert abs(rc.correlation() - _batch_corr(xs, ys)) < 1e-9
        assert abs(rc.variance_x() - statistics.variance(xs)) < 1e-7


def test_perfect_correlation():
    xs = list(range(100))
    assert abs(RunningCovariance(xs, [2 * x + 5 for x in xs]).correlation() - 1.0) < 1e-9
    assert abs(RunningCovariance(xs, [-3 * x for x in xs]).correlation() + 1.0) < 1e-9


def test_merge_exactness():
    rng = random.Random(2)
    for _ in range(300):
        n = rng.randint(2, 40)
        xs = [rng.uniform(-5, 5) for _ in range(n)]
        ys = [rng.uniform(-5, 5) for _ in range(n)]
        k = rng.randint(1, n - 1)
        merged = RunningCovariance(xs[:k], ys[:k]) + RunningCovariance(xs[k:], ys[k:])
        full = RunningCovariance(xs, ys)
        assert merged.n == n
        assert abs(merged.covariance() - full.covariance()) < 1e-7
        assert abs(merged.correlation() - full.correlation()) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        RunningCovariance([1], [1]).covariance()
    with pytest.raises(ValueError):
        RunningCovariance([1, 1], [5, 5]).correlation()
    with pytest.raises(ValueError):
        RunningCovariance([1, 2], [1])
