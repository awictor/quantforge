"""FFT cross-correlation and lead-lag detection."""

import math

import pytest

from quantforge import (
    cross_correlation,
    normalized_cross_correlation,
    lag_at_max_correlation,
)


def _brute(x, y, lag):
    n = len(x)
    return sum(x[k + lag] * y[k] for k in range(n) if 0 <= k + lag < n)


def test_fft_matches_brute_force():
    N = 64
    x = [math.sin(2 * math.pi * 3 * n / N) + 0.3 * n / N for n in range(N)]
    y = [math.cos(2 * math.pi * 2 * n / N) for n in range(N)]
    lags, vals = cross_correlation(x, y, max_lag=20)
    assert max(abs(vals[i] - _brute(x, y, lags[i])) for i in range(len(lags))) < 1e-10


def test_recovers_known_lag():
    N = 64
    d = 5
    base = [math.sin(2 * math.pi * 4 * n / N) for n in range(N)]
    x = [base[(n - d) % N] for n in range(N)]   # x is base delayed by d
    lag, coef = lag_at_max_correlation(x, base, max_lag=20)
    assert lag == d
    assert coef > 0.8


def test_identical_series_peaks_at_zero():
    base = [math.sin(2 * math.pi * 4 * n / 64) for n in range(64)]
    lag, coef = lag_at_max_correlation(base, base, max_lag=10)
    assert lag == 0
    assert abs(coef - 1.0) < 1e-9


def test_normalized_in_range():
    N = 64
    x = [math.sin(2 * math.pi * 3 * n / N) for n in range(N)]
    y = [math.cos(2 * math.pi * 2 * n / N) for n in range(N)]
    _, nv = normalized_cross_correlation(x, y, max_lag=20)
    assert all(-1.0001 <= v <= 1.0001 for v in nv)


def test_lags_are_symmetric():
    x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    lags, vals = cross_correlation(x, x, max_lag=3)
    assert lags == [-3, -2, -1, 0, 1, 2, 3]
    # autocorrelation is symmetric about lag 0
    assert abs(vals[2] - vals[4]) < 1e-9
    assert abs(vals[1] - vals[5]) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        cross_correlation([1.0, 2.0], [1.0])         # length mismatch
    with pytest.raises(ValueError):
        cross_correlation([], [])                     # empty
    with pytest.raises(ValueError):
        cross_correlation([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0], max_lag=4)
    with pytest.raises(ValueError):
        normalized_cross_correlation([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])  # zero variance
