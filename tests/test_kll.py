import random

import pytest

from quantforge import KLL


def test_small_stream_quantiles():
    # fast coverage without the 200k-sample sketches
    random.seed(3)
    data = [random.gauss(0, 1) for _ in range(20000)]
    kll = KLL(k=200, seed=1).add_all(data)
    s = sorted(data)
    assert kll.n == len(data)
    assert abs(kll.quantile(0.5) - s[len(s) // 2]) < 0.1


def _exact_rank(data, x):
    return sum(1 for v in data if v <= x)


def _data(n=200000, seed=0):
    random.seed(seed)
    return [random.gauss(0, 1) for _ in range(n)]


@pytest.mark.slow
def test_rank_error_is_bounded():
    data = _data()
    kll = KLL(k=200, seed=1).add_all(data)
    assert kll.n == len(data)
    maxerr = max(abs(kll.rank(x) - _exact_rank(data, x)) / len(data)
                 for x in (-2, -1, -0.5, 0, 0.5, 1, 2, 3))
    assert maxerr < 0.02


@pytest.mark.slow
def test_bulk_quantiles_close_in_value():
    data = _data()
    kll = KLL(k=200, seed=1).add_all(data)
    s = sorted(data)
    for q in (0.1, 0.25, 0.5, 0.75, 0.9):
        ex = s[min(len(s) - 1, int(q * len(s)))]
        assert abs(kll.quantile(q) - ex) < 0.08


@pytest.mark.slow
def test_tail_quantiles_rank_consistent():
    data = _data()
    kll = KLL(k=200, seed=1).add_all(data)
    for q in (0.01, 0.99):
        v = kll.quantile(q)
        assert abs(_exact_rank(data, v) - q * len(data)) / len(data) < 0.02


@pytest.mark.slow
def test_cdf():
    kll = KLL(k=200, seed=1).add_all(_data())
    assert kll.cdf(-10) < 0.01
    assert kll.cdf(10) > 0.99
    assert abs(kll.cdf(0) - 0.5) < 0.02


def test_reproducible():
    data = _data(50000)
    a = KLL(200, seed=7).add_all(data)
    b = KLL(200, seed=7).add_all(data)
    assert a.quantile(0.5) == b.quantile(0.5)


@pytest.mark.slow
def test_merge():
    random.seed(1)
    d1 = [random.gauss(-1, 1) for _ in range(100000)]
    d2 = [random.gauss(2, 1) for _ in range(100000)]
    m = KLL(200, seed=3).add_all(d1)
    m.merge(KLL(200, seed=4).add_all(d2))
    comb = sorted(d1 + d2)
    assert m.n == 200000
    for q in (0.25, 0.5, 0.75):
        assert abs(m.quantile(q) - comb[int(q * len(comb))]) < 0.15


@pytest.mark.slow
def test_compact_memory():
    kll = KLL(k=200, seed=1).add_all(_data())
    assert sum(len(l) for l in kll._levels) < 2000


def test_empty_raises():
    with pytest.raises(ValueError):
        KLL().quantile(0.5)
