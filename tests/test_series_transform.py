"""Series acceleration: Wynn's epsilon and Euler transform."""

import math

import pytest

from quantforge import wynn_epsilon, euler_transform


def _leibniz_ps(n):
    s, acc = [], 0.0
    for k in range(n):
        acc += (-1) ** k / (2 * k + 1)
        s.append(acc)
    return s


def test_wynn_leibniz_machine_precision():
    w = wynn_epsilon(_leibniz_ps(20))
    assert abs(w * 4 - math.pi) < 1e-12


def test_wynn_exact_on_geometric():
    s, acc = [], 0.0
    for k in range(8):
        acc += 0.5 ** k
        s.append(acc)
    assert abs(wynn_epsilon(s) - 2.0) < 1e-12


def test_wynn_beats_raw_on_zeta2():
    s, acc = [], 0.0
    for k in range(1, 30):
        acc += 1 / k ** 2
        s.append(acc)
    w = wynn_epsilon(s)
    raw_err = abs(s[-1] - math.pi ** 2 / 6)
    assert abs(w - math.pi ** 2 / 6) < raw_err / 2


def test_euler_leibniz():
    e = euler_transform([1 / (2 * k + 1) for k in range(25)])
    assert abs(e * 4 - math.pi) < 1e-6


def test_euler_ln2():
    e = euler_transform([1 / (k + 1) for k in range(30)])
    assert abs(e - math.log(2)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        wynn_epsilon([1.0, 2.0])
    with pytest.raises(ValueError):
        euler_transform([])
