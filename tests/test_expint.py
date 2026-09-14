"""Tests for exponential integrals E1, Ei, E_n, cross-checked against numerical integration."""

import math

import pytest

from quantforge.expint import e1, ei, en


def _close(a, b, tol=1e-8):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _num_e1(x, N=200000):
    a, b = x, x + 50
    h = (b - a) / N
    total = 0.0
    for i in range(N + 1):
        t = a + i * h
        w = 0.5 if i in (0, N) else 1.0
        total += w * math.exp(-t) / t
    return total * h


def _num_en(n, x, N=100000):
    a, b = 1.0, 60.0
    h = (b - a) / N
    total = 0.0
    for i in range(N + 1):
        t = a + i * h
        w = 0.5 if i in (0, N) else 1.0
        total += w * math.exp(-x * t) / t ** n
    return total * h


def test_e1_known_values():
    assert _close(e1(1.0), 0.21938393439552027)
    assert _close(e1(2.0), 0.04890051070806112)
    assert _close(e1(0.5), 0.5597735947761608)


def test_ei_known_values():
    assert _close(ei(1.0), 1.8951178163559366)
    assert _close(ei(2.0), 4.954234356001890)
    assert _close(ei(0.5), 0.45421990486584387)


def test_e1_vs_integration():
    for x in (0.3, 0.7, 1.5, 3.0, 8.0):
        assert _close(e1(x), _num_e1(x), 1e-6)


def test_en_vs_integration():
    for n in (2, 3, 5):
        for x in (0.5, 1.0, 3.0):
            assert _close(en(n, x), _num_en(n, x), 1e-5)


def test_en_at_zero():
    assert _close(en(2, 0), 1.0)
    assert _close(en(3, 0), 0.5)
    assert _close(en(5, 0), 0.25)


def test_e1_equals_en_order_one():
    for x in (0.3, 1.0, 5.0):
        assert _close(en(1, x), e1(x))


def test_ei_negative_relation():
    for x in (0.5, 2.0, 6.0):
        assert _close(ei(-x), -e1(x))


def test_ei_large_x_asymptotic():
    # continuity across the series/asymptotic boundary and against exp(x)/x leading term
    x = 50.0
    assert ei(x) > 0
    # leading asymptotic term is exp(x)/x; the full value adds ~1/x corrections
    assert _close(ei(x) * x / math.exp(x), 1.0, 0.05)


def test_e1_monotone_decreasing():
    xs = [0.5, 1.0, 2.0, 5.0, 10.0]
    vals = [e1(x) for x in xs]
    assert all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))


def test_e0_closed_form():
    for x in (0.5, 1.0, 3.0):
        assert _close(en(0, x), math.exp(-x) / x)


def test_domain_errors():
    with pytest.raises(ValueError):
        e1(0)
    with pytest.raises(ValueError):
        e1(-1)
    with pytest.raises(ValueError):
        ei(0)
    with pytest.raises(ValueError):
        en(-1, 1)
    with pytest.raises(ValueError):
        en(0, 0)
