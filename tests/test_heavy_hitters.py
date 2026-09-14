"""Tests for streaming heavy-hitters: MisraGries and SpaceSaving.

Each estimate is cross-checked against the exact Counter over the same stream, verifying
the algorithms' guarantees rather than exact frequencies.
"""

import random
from collections import Counter

import pytest

from quantforge.heavy_hitters import MisraGries, SpaceSaving


def _skewed_stream(seed=1, n=20000):
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        r = rng.random()
        if r < 0.30:
            out.append("A")
        elif r < 0.50:
            out.append("B")
        elif r < 0.62:
            out.append("C")
        else:
            out.append("t%d" % rng.randint(0, 400))
    rng.shuffle(out)
    return out


def test_misra_gries_never_overcounts():
    stream = _skewed_stream()
    exact = Counter(stream)
    mg = MisraGries(10)
    for x in stream:
        mg.add(x)
    for it, c in mg.counts().items():
        assert c <= exact[it]  # stored count never exceeds the true count


def test_misra_gries_undercount_bounded():
    stream = _skewed_stream()
    exact = Counter(stream)
    k = 10
    mg = MisraGries(k)
    for x in stream:
        mg.add(x)
    bound = len(stream) / (k + 1)
    for it, c in mg.counts().items():
        assert exact[it] - c <= bound + 1e-9


def test_misra_gries_keeps_true_heavy_hitters():
    stream = _skewed_stream()
    exact = Counter(stream)
    k = 10
    mg = MisraGries(k)
    for x in stream:
        mg.add(x)
    bound = len(stream) / (k + 1)
    for it, c in exact.items():
        if c > bound:  # guaranteed to survive
            assert it in mg.counts()


def test_misra_gries_counter_cap():
    stream = _skewed_stream()
    mg = MisraGries(10)
    for x in stream:
        mg.add(x)
    assert len(mg.counts()) <= 10


def test_misra_gries_exact_when_distinct_fits():
    mg = MisraGries(5)
    for x in "aaabbc":
        mg.add(x)
    # only 3 distinct items, all fit -> exact counts
    assert mg.counts() == {"a": 3, "b": 2, "c": 1}


def test_misra_gries_weighted_matches_repeated():
    a = MisraGries(5)
    a.add("x", 4)
    a.add("y", 2)
    b = MisraGries(5)
    for _ in range(4):
        b.add("x")
    for _ in range(2):
        b.add("y")
    assert a.counts() == b.counts()
    assert a.n == b.n == 6


def test_misra_gries_heavy_hitters_threshold():
    mg = MisraGries(5)
    for x in "aaaaaaabbbc":  # a=7, b=3, c=1, n=11
        mg.add(x)
    hh = mg.heavy_hitters(0.5)  # >= 5.5
    assert set(hh) == {"a"}


def test_space_saving_count_is_upper_bound():
    stream = _skewed_stream()
    exact = Counter(stream)
    ss = SpaceSaving(10)
    for x in stream:
        ss.add(x)
    counts = ss.counts()
    for it, c in counts.items():
        assert c >= exact[it]  # Space-Saving overcounts, never under


def test_space_saving_lower_bound_holds():
    stream = _skewed_stream()
    exact = Counter(stream)
    ss = SpaceSaving(10)
    for x in stream:
        ss.add(x)
    counts = ss.counts()
    errors = ss.errors()
    for it in counts:
        # count - error <= true <= count
        assert counts[it] - errors[it] <= exact[it] <= counts[it]


def test_space_saving_top_k_matches_exact_on_skew():
    stream = _skewed_stream()
    exact = Counter(stream)
    ss = SpaceSaving(10)
    for x in stream:
        ss.add(x)
    ss_top = [it for it, _ in ss.top(3)]
    exact_top = [it for it, _ in exact.most_common(3)]
    assert ss_top == exact_top


def test_space_saving_guaranteed_lower_bounds_positive_for_heavy():
    stream = _skewed_stream()
    ss = SpaceSaving(10)
    for x in stream:
        ss.add(x)
    g = dict(ss.guaranteed(3))
    # the three dominant items have large guaranteed lower bounds
    for it in ("A", "B", "C"):
        assert g[it] > 0


def test_space_saving_exact_when_distinct_fits():
    ss = SpaceSaving(5)
    for x in "aaabbc":
        ss.add(x)
    assert ss.counts() == {"a": 3, "b": 2, "c": 1}
    assert all(e == 0 for e in ss.errors().values())


def test_space_saving_eviction_sets_error():
    ss = SpaceSaving(2)
    # fill with a, b; then c evicts the min and inherits its count as error
    for x in "aabc":  # a=2, b=1, then c
        ss.add(x)
    counts = ss.counts()
    errors = ss.errors()
    assert "a" in counts  # a is heaviest, stays
    assert len(counts) == 2
    # the newest survivor carries a nonzero error from the eviction
    assert max(errors.values()) >= 1


def test_counter_cap_space_saving():
    stream = _skewed_stream()
    ss = SpaceSaving(8)
    for x in stream:
        ss.add(x)
    assert len(ss.counts()) <= 8


@pytest.mark.parametrize("cls", [MisraGries, SpaceSaving])
def test_k_below_one_raises(cls):
    with pytest.raises(ValueError):
        cls(0)


@pytest.mark.parametrize("cls", [MisraGries, SpaceSaving])
def test_nonpositive_weight_raises(cls):
    with pytest.raises(ValueError):
        cls(3).add("x", 0)


def test_misra_gries_threshold_out_of_range_raises():
    with pytest.raises(ValueError):
        MisraGries(3).heavy_hitters(1.5)


def test_space_saving_negative_top_raises():
    with pytest.raises(ValueError):
        SpaceSaving(3).top(-1)


def test_add_returns_self():
    mg = MisraGries(3)
    assert mg.add("x") is mg
    ss = SpaceSaving(3)
    assert ss.add("x") is ss
