"""Streaming mean/variance/skewness/kurtosis (Terriberry)."""

import random
import statistics

import pytest

from quantforge import RunningMoments


def _pop_skew(d):
    n = len(d); m = sum(d) / n
    m2 = sum((x - m) ** 2 for x in d) / n
    m3 = sum((x - m) ** 3 for x in d) / n
    return m3 / m2 ** 1.5


def _pop_kurt(d):
    n = len(d); m = sum(d) / n
    m2 = sum((x - m) ** 2 for x in d) / n
    m4 = sum((x - m) ** 4 for x in d) / n
    return m4 / m2 ** 2 - 3


def test_matches_batch_statistics():
    rng = random.Random(3)
    data = [rng.gauss(5, 2) for _ in range(2000)]
    rm = RunningMoments(data)
    assert abs(rm.mean - statistics.mean(data)) < 1e-9
    assert abs(rm.variance() - statistics.variance(data)) < 1e-9
    assert abs(rm.skewness() - _pop_skew(data)) < 1e-6
    assert abs(rm.kurtosis() - _pop_kurt(data)) < 1e-6


def test_merge_equals_full():
    rng = random.Random(5)
    data = [rng.gauss(0, 1) for _ in range(1000)]
    a = RunningMoments(data[:700])
    b = RunningMoments(data[700:])
    merged = a + b
    full = RunningMoments(data)
    assert abs(merged.mean - full.mean) < 1e-9
    assert abs(merged.M2 - full.M2) < 1e-6
    assert abs(merged.M3 - full.M3) < 1e-6
    assert abs(merged.M4 - full.M4) < 1e-6


def test_multiway_merge_associative():
    rng = random.Random(7)
    data = [rng.gauss(2, 3) for _ in range(2000)]
    full = RunningMoments(data)
    chunks = [RunningMoments(data[i * 200:(i + 1) * 200]) for i in range(10)]
    acc = chunks[0]
    for c in chunks[1:]:
        acc = acc + c
    assert abs(acc.M4 - full.M4) < 1e-6
    assert abs(acc.skewness() - full.skewness()) < 1e-9


def test_symmetric_zero_skew():
    assert abs(RunningMoments([-3, -2, -1, 0, 1, 2, 3]).skewness()) < 1e-12


def test_normal_zero_excess_kurtosis():
    rng = random.Random(9)
    rm = RunningMoments([rng.gauss(0, 1) for _ in range(50000)])
    assert abs(rm.kurtosis()) < 0.1


def test_merge_with_empty():
    a = RunningMoments([1, 2, 3])
    empty = RunningMoments()
    assert (a + empty).mean == a.mean
    assert (empty + a).M2 == a.M2
