"""Tests for Airy functions Ai, Bi, cross-checked against known values, the ODE, and Wronskian."""

import math

import pytest

from quantforge.airy import airy_ai, airy_bi


def _close(a, b, tol=1e-8):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _deriv(fn, x, h=1e-5):
    return (fn(x + h) - fn(x - h)) / (2 * h)


def test_known_values():
    assert _close(airy_ai(0), 0.3550280538878172)
    assert _close(airy_bi(0), 0.6149266274460007)
    assert _close(airy_ai(1), 0.13529241631288147, 1e-7)
    assert _close(airy_bi(1), 1.2074235949528713, 1e-7)
    assert _close(airy_ai(-1), 0.5355608832923521, 1e-7)
    assert _close(airy_bi(-1), 0.10399738949694459, 1e-6)
    assert _close(airy_ai(2), 0.03492413042327235, 1e-7)


def test_ode_y_double_prime_equals_x_y():
    for fn in (airy_ai, airy_bi):
        for x in (-2, -0.5, 0.5, 1.5, 3):
            h = 1e-4
            ypp = (fn(x + h) - 2 * fn(x) + fn(x - h)) / h ** 2
            assert _close(ypp, x * fn(x), 1e-3)


def test_wronskian():
    for x in (-1, 0, 1, 2):
        w = airy_ai(x) * _deriv(airy_bi, x) - _deriv(airy_ai, x) * airy_bi(x)
        assert _close(w, 1 / math.pi, 1e-4)


def test_first_zero_of_ai():
    assert abs(airy_ai(-2.33810741045976)) < 1e-6


def test_decay_and_growth():
    # Ai decays toward 0 for large x; Bi grows
    assert 0 < airy_ai(5) < airy_ai(2)
    assert airy_bi(5) > airy_bi(2) > airy_bi(0)


def test_negative_oscillation():
    # both oscillate for x < 0
    a = airy_ai(-5)
    assert -1 < a < 1


def test_ai_positive_small():
    assert airy_ai(0.5) > 0
    assert airy_bi(0.5) > 0


def test_second_zero_of_ai():
    assert abs(airy_ai(-4.08794944413097)) < 1e-6


def test_ratio_at_zero():
    # Bi(0)/Ai(0) = sqrt(3)
    assert _close(airy_bi(0) / airy_ai(0), math.sqrt(3), 1e-9)


def test_moderate_negative_bounded():
    # amplitude ~ |x|^{-1/4} / sqrt(pi); stays small in the series' reliable range
    for x in (-5, -8, -10):
        assert abs(airy_ai(x)) < 1.0
        assert abs(airy_bi(x)) < 1.0
