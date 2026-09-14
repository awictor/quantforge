"""Find all roots of a function on an interval by sign-change scanning."""

import math

import pytest

from quantforge import find_all_roots, count_sign_changes


def test_cubic_roots():
    f = lambda x: (x - 1) * (x - 2) * (x - 3)
    r = find_all_roots(f, 0, 4)
    assert len(r) == 3
    assert all(abs(f(v)) < 1e-8 for v in r)
    for got, want in zip(r, [1, 2, 3]):
        assert abs(got - want) < 1e-6


def test_sine_roots():
    r = find_all_roots(math.sin, 0.5, 10)
    assert len(r) == 3
    for k, v in enumerate(r, 1):
        assert abs(v - k * math.pi) < 1e-6


def test_count_sign_changes():
    assert count_sign_changes(math.sin, 0.5, 10) == 3
    assert count_sign_changes(math.cos, 0, 20) == 6


def test_no_roots_and_node_zero():
    assert find_all_roots(lambda x: x * x + 1, -5, 5) == []
    assert find_all_roots(lambda x: x, -1, 1) == [0.0]


def test_transcendental_roots():
    f = lambda x: math.exp(x) - 3 * x
    r = find_all_roots(f, 0, 3)
    assert len(r) == 2
    assert max(abs(f(v)) for v in r) < 1e-8


def test_validation():
    with pytest.raises(ValueError):
        find_all_roots(math.sin, 5, 5)
    with pytest.raises(ValueError):
        count_sign_changes(math.sin, 0, 1, n=0)
