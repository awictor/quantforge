import math

import pytest

from quantforge import bdf2, backward_euler


def close(a, b, tol=1e-2):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _exact(t):
    return (t / 2 - 0.25) + 1.25 * math.exp(-2 * t)


def _f(t, y):
    return -2 * y + t


def test_decay():
    _, ys = bdf2(lambda t, y: -y, 1.0, 0, 2, 100)
    assert close(ys[-1][0], math.exp(-2), 1e-3)


def test_second_order_convergence():
    _, y1 = bdf2(_f, 1.0, 0, 1, 20)
    _, y2 = bdf2(_f, 1.0, 0, 1, 40)
    e1 = abs(y1[-1][0] - _exact(1))
    e2 = abs(y2[-1][0] - _exact(1))
    assert 3.0 < e1 / e2 < 5.0


def test_beats_backward_euler():
    _, yb2 = bdf2(_f, 1.0, 0, 1, 40)
    _, ybe = backward_euler(_f, 1.0, 0, 1, 40)
    assert abs(yb2[-1][0] - _exact(1)) < abs(ybe[-1][0] - _exact(1))


def test_stiff_stability():
    def stiff(t, y):
        return -1000 * (y - math.cos(t)) - math.sin(t)

    _, ys = bdf2(stiff, 0.0, 0, 3, 60)
    assert abs(ys[-1][0] - math.cos(3)) < 0.05


def test_vector_harmonic():
    def harm(t, y):
        return [y[1], -y[0]]

    _, ys = bdf2(harm, [1.0, 0.0], 0, math.pi, 400)
    assert close(ys[-1][0], -1.0, 2e-2)
    assert abs(ys[-1][1]) < 2e-2


def test_too_few_steps_raises():
    with pytest.raises(ValueError):
        bdf2(_f, 1.0, 0, 1, 1)
