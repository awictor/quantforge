import math

from quantforge import backward_euler, trapezoidal


def close(a, b, tol=1e-2):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_backward_euler_decay():
    _, ys = backward_euler(lambda t, y: -y, 1.0, 0, 2, 200)
    assert close(ys[-1][0], math.exp(-2), 0.02)


def test_trapezoidal_decay_more_accurate():
    _, ys = trapezoidal(lambda t, y: -y, 1.0, 0, 2, 50)
    assert close(ys[-1][0], math.exp(-2), 1e-3)


def test_stiff_stability():
    # y' = -1000(y - cos t) - sin t : explicit RK needs h < 0.002, implicit is stable at h=0.05
    def stiff(t, y):
        return -1000 * (y - math.cos(t)) - math.sin(t)

    _, ys = backward_euler(stiff, 0.0, 0, 3, 60)
    assert abs(ys[-1][0] - math.cos(3)) < 0.05


def _exact(t):
    return (t / 2 - 0.25) + 1.25 * math.exp(-2 * t)


def _f(t, y):
    return -2 * y + t


def test_trapezoidal_second_order():
    _, y1 = trapezoidal(_f, 1.0, 0, 1, 20)
    _, y2 = trapezoidal(_f, 1.0, 0, 1, 40)
    e1 = abs(y1[-1][0] - _exact(1))
    e2 = abs(y2[-1][0] - _exact(1))
    assert 3.0 < e1 / e2 < 5.0     # ~4 for order 2


def test_backward_euler_first_order():
    _, y1 = backward_euler(_f, 1.0, 0, 1, 20)
    _, y2 = backward_euler(_f, 1.0, 0, 1, 40)
    e1 = abs(y1[-1][0] - _exact(1))
    e2 = abs(y2[-1][0] - _exact(1))
    assert 1.7 < e1 / e2 < 2.3     # ~2 for order 1


def test_vector_harmonic_oscillator():
    def harm(t, y):
        return [y[1], -y[0]]

    _, ys = trapezoidal(harm, [1.0, 0.0], 0, math.pi, 200)
    assert close(ys[-1][0], -1.0, 1e-2)
    assert abs(ys[-1][1]) < 1e-2
