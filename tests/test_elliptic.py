"""Tests for complete elliptic integrals and AGM, cross-checked against integration."""

import math
import random

import pytest

from quantforge.elliptic import agm, elliptic_k, elliptic_e


def _close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _num_k(m, N=200000):
    h = (math.pi / 2) / N
    s = 0.0
    for i in range(N + 1):
        th = i * h
        w = 0.5 if i in (0, N) else 1.0
        s += w / math.sqrt(1 - m * math.sin(th) ** 2)
    return s * h


def _num_e(m, N=200000):
    h = (math.pi / 2) / N
    s = 0.0
    for i in range(N + 1):
        th = i * h
        w = 0.5 if i in (0, N) else 1.0
        s += w * math.sqrt(1 - m * math.sin(th) ** 2)
    return s * h


def test_known_values():
    assert _close(elliptic_k(0), math.pi / 2)
    assert _close(elliptic_e(0), math.pi / 2)
    assert _close(elliptic_e(1), 1.0)
    assert _close(elliptic_k(0.5), _num_k(0.5), 1e-9)


def test_vs_integration_fast():
    for m in (0.1, 0.4, 0.75, 0.95):
        assert _close(elliptic_k(m), _num_k(m, N=40000), 1e-5)
        assert _close(elliptic_e(m), _num_e(m, N=40000), 1e-5)


@pytest.mark.slow
def test_fuzz_vs_integration():
    rng = random.Random(551)
    for _ in range(500):
        m = rng.uniform(0, 0.98)
        assert _close(elliptic_k(m), _num_k(m), 1e-6)
        assert _close(elliptic_e(m), _num_e(m), 1e-6)


def test_legendre_relation():
    for m in (0.2, 0.5, 0.7, 0.9):
        lhs = (elliptic_e(m) * elliptic_k(1 - m) + elliptic_e(1 - m) * elliptic_k(m)
               - elliptic_k(m) * elliptic_k(1 - m))
        assert _close(lhs, math.pi / 2, 1e-9)


def test_agm_properties():
    assert _close(agm(1, 1), 1.0)
    assert _close(agm(5, 5), 5.0)
    assert _close(agm(2, 8), agm(8, 2))  # symmetric
    assert _close(agm(24, 6), 13.4581714817256, 1e-9)


def test_agm_between_gm_and_am():
    rng = random.Random(552)
    for _ in range(1000):
        a = rng.uniform(0.1, 100)
        b = rng.uniform(0.1, 100)
        gm = math.sqrt(a * b)
        am = (a + b) / 2
        v = agm(a, b)
        assert min(gm, am) - 1e-9 <= v <= max(gm, am) + 1e-9


def test_k_increases_with_m():
    ms = [0.0, 0.3, 0.6, 0.9, 0.99]
    ks = [elliptic_k(m) for m in ms]
    assert all(ks[i] < ks[i + 1] for i in range(len(ks) - 1))


def test_e_decreases_with_m():
    ms = [0.0, 0.3, 0.6, 0.9, 1.0]
    es = [elliptic_e(m) for m in ms]
    assert all(es[i] > es[i + 1] for i in range(len(es) - 1))


def test_domain_errors():
    with pytest.raises(ValueError):
        elliptic_k(1)
    with pytest.raises(ValueError):
        elliptic_k(-0.1)
    with pytest.raises(ValueError):
        elliptic_e(1.5)
    with pytest.raises(ValueError):
        agm(-1, 2)


def test_pendulum_period_sanity():
    # small-amplitude pendulum period ratio T/T0 = (2/pi) K(sin^2(theta0/2))
    theta0 = 0.1
    m = math.sin(theta0 / 2) ** 2
    ratio = (2 / math.pi) * elliptic_k(m)
    assert ratio > 1.0  # always longer than the linear period
    assert _close(ratio, 1.0, 1e-2)  # but only slightly, for small amplitude
