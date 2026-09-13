"""Sequence acceleration: Aitken, Shanks, Steffensen."""

import math

import pytest

from quantforge import aitken, shanks, steffensen


def _leibniz(n):
    s, acc = [], 0.0
    for k in range(n):
        acc += (-1) ** k / (2 * k + 1)
        s.append(acc)
    return s


def test_aitken_speeds_up_leibniz():
    seq = _leibniz(20)
    acc = aitken(seq)
    raw_err = abs(seq[-1] * 4 - math.pi)
    acc_err = abs(acc[-1] * 4 - math.pi)
    assert acc_err < raw_err / 100


def test_aitken_exact_on_geometric():
    r = 0.5
    s, acc = [], 0.0
    for k in range(6):
        acc += r ** k
        s.append(acc)
    a = aitken(s)
    assert abs(a[0] - 2.0) < 1e-12         # sum r^k = 1/(1-r) = 2


def test_shanks_equals_aitken():
    seq = _leibniz(15)
    assert shanks(seq) == aitken(seq)


def test_steffensen_cos_fixed_point():
    res = steffensen(math.cos, 0.5)
    assert res["converged"]
    assert abs(math.cos(res["root"]) - res["root"]) < 1e-10
    assert res["iterations"] < 15          # far fewer than plain iteration (~69)


def test_steffensen_sqrt2():
    res = steffensen(lambda x: (x + 2 / x) / 2, 1.0)
    assert abs(res["root"] - math.sqrt(2)) < 1e-12


def test_steffensen_beats_plain_iteration():
    def plain(g, x0, tol=1e-12, maxit=1000):
        x = x0
        for i in range(1, maxit + 1):
            xn = g(x)
            if abs(xn - x) < tol:
                return i
            x = xn
        return maxit
    assert steffensen(math.cos, 0.5)["iterations"] < plain(math.cos, 0.5) / 4


def test_validation():
    with pytest.raises(ValueError):
        aitken([1.0, 2.0])
