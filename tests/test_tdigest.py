import random

import pytest

from quantforge import TDigest


def _exact_q(data, q):
    s = sorted(data)
    if q <= 0:
        return s[0]
    if q >= 1:
        return s[-1]
    idx = q * (len(s) - 1)
    lo = int(idx)
    frac = idx - lo
    return s[lo] if lo + 1 >= len(s) else s[lo] + frac * (s[lo + 1] - s[lo])


def test_quantile_accuracy_gaussian():
    random.seed(0)
    data = [random.gauss(0, 1) for _ in range(100000)]
    td = TDigest(compression=200).add_all(data)
    assert abs(td.count - 100000) < 1e-6
    for q in (0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999):
        approx = td.quantile(q)
        ex = _exact_q(data, q)
        assert abs(approx - ex) < 0.07


def test_median_and_p99():
    random.seed(0)
    td = TDigest(200).add_all([random.gauss(0, 1) for _ in range(100000)])
    assert abs(td.quantile(0.5)) < 0.02
    assert abs(td.quantile(0.99) - 2.326) < 0.06


def test_cdf_is_inverse_and_monotone():
    random.seed(0)
    td = TDigest(200).add_all([random.gauss(0, 1) for _ in range(100000)])
    assert abs(td.cdf(0.0) - 0.5) < 0.02
    assert td.cdf(-10) < 0.01
    assert td.cdf(10) > 0.99


def test_merge_equals_combined_stream():
    random.seed(1)
    d1 = [random.gauss(-1, 1) for _ in range(50000)]
    d2 = [random.gauss(3, 2) for _ in range(50000)]
    a = TDigest(200).add_all(d1)
    a.merge(TDigest(200).add_all(d2))
    combined = TDigest(200).add_all(d1 + d2)
    for q in (0.1, 0.5, 0.9, 0.99):
        assert abs(a.quantile(q) - combined.quantile(q)) < 0.1
    assert abs(a.count - 100000) < 1e-6


def test_uniform_quantiles():
    random.seed(2)
    td = TDigest(200).add_all([random.random() for _ in range(100000)])
    for q in (0.1, 0.3, 0.5, 0.7, 0.9):
        assert abs(td.quantile(q) - q) < 0.02


def test_empty_raises():
    with pytest.raises(ValueError):
        TDigest().quantile(0.5)
