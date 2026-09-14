"""Tests for the Lambert W function, cross-checked via the defining identity w e^w = x."""

import math
import random

import pytest

from quantforge.lambert import lambert_w0, lambert_wm1

_INV_E = 1.0 / math.e


def _close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def test_known_values():
    assert _close(lambert_w0(0), 0.0)
    assert _close(lambert_w0(math.e), 1.0)
    assert _close(lambert_w0(1), 0.5671432904097838)  # the omega constant
    assert _close(lambert_w0(-_INV_E), -1.0, 1e-6)
    assert _close(lambert_wm1(-_INV_E), -1.0, 1e-6)


def test_fuzz_w0_identity():
    rng = random.Random(541)
    for _ in range(10000):
        if rng.random() < 0.7:
            x = rng.uniform(-_INV_E + 1e-9, 1000)
        else:
            x = rng.uniform(-_INV_E + 1e-9, 0)
        w = lambert_w0(x)
        assert _close(w * math.exp(w), x, 1e-7)
        assert w >= -1 - 1e-9


def test_fuzz_wm1_identity():
    rng = random.Random(542)
    for _ in range(10000):
        x = rng.uniform(-_INV_E + 1e-9, -1e-6)
        w = lambert_wm1(x)
        assert _close(w * math.exp(w), x, 1e-6)
        assert w <= -1 + 1e-6


def test_large_x():
    for x in (1e3, 1e6, 1e9):
        w = lambert_w0(x)
        assert _close(w * math.exp(w), x, 1e-5)


def test_positive_branch_monotone():
    xs = [0.1, 1.0, 10.0, 100.0]
    vals = [lambert_w0(x) for x in xs]
    assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))


def test_negative_domain_both_branches_distinct():
    x = -0.2
    w0 = lambert_w0(x)
    wm1 = lambert_wm1(x)
    assert w0 >= -1
    assert wm1 <= -1
    assert _close(w0 * math.exp(w0), x, 1e-7)
    assert _close(wm1 * math.exp(wm1), x, 1e-6)


def test_solves_exponential_equation():
    # solve x e^x = 5 -> x = W0(5)
    w = lambert_w0(5.0)
    assert _close(w * math.exp(w), 5.0, 1e-9)


def test_w0_of_small_negative():
    x = -0.001
    w = lambert_w0(x)
    assert _close(w * math.exp(w), x, 1e-9)
    assert w < 0


def test_w0_below_branch_point_raises():
    with pytest.raises(ValueError):
        lambert_w0(-1.0)


def test_wm1_out_of_domain_raises():
    with pytest.raises(ValueError):
        lambert_wm1(1.0)
    with pytest.raises(ValueError):
        lambert_wm1(-1.0)  # below -1/e
