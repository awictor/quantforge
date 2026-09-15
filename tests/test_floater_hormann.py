import math

import pytest

from quantforge import floater_hormann_weights, floater_hormann_interpolate as fh
from quantforge.barycentric import barycentric_weights, barycentric_eval


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_interpolates_nodes_exactly():
    xs = [i * 0.5 for i in range(9)]
    ys = [math.sin(x) for x in xs]
    for d in (0, 3, 5):
        for k in range(len(xs)):
            assert close(fh(xs, ys, xs[k], d), ys[k])


def test_exact_for_low_degree_polynomial():
    poly = lambda x: 2 * x ** 3 - x ** 2 + 3 * x - 1
    xs = [i * 0.3 for i in range(8)]
    ys = [poly(x) for x in xs]
    for xq in (0.15, 0.9, 1.7):
        assert close(fh(xs, ys, xq, d=3), poly(xq), 1e-6)


def test_runge_controlled_vs_polynomial():
    runge = lambda x: 1 / (1 + 25 * x * x)
    N = 21
    xs = [-1 + 2 * i / (N - 1) for i in range(N)]
    ys = [runge(x) for x in xs]
    fh_max = max(abs(fh(xs, ys, -1 + 2 * i / 2000, d=3) - runge(-1 + 2 * i / 2000))
                 for i in range(2001))
    pw = barycentric_weights(xs)
    poly_max = max(abs(barycentric_eval(xs, ys, pw, -1 + 2 * i / 2000) - runge(-1 + 2 * i / 2000))
                   for i in range(2001))
    assert fh_max < poly_max
    assert fh_max < 0.1


def test_no_poles_on_grid():
    runge = lambda x: 1 / (1 + 25 * x * x)
    xs = [-1 + 2 * i / 20 for i in range(21)]
    ys = [runge(x) for x in xs]
    assert all(abs(fh(xs, ys, -1 + 2 * i / 5000, d=3)) < 1e6 for i in range(5001))


def test_berrut_weights_alternate():
    w = floater_hormann_weights([0, 1, 2, 3, 4], d=0)
    signs = [1 if x > 0 else -1 for x in w]
    assert all(signs[i] != signs[i + 1] for i in range(len(signs) - 1))


def test_higher_d_improves_accuracy():
    xs = [-1 + 2 * i / 20 for i in range(21)]
    ys = [math.exp(x) for x in xs]
    errs = {d: max(abs(fh(xs, ys, -1 + 2 * i / 2000, d) - math.exp(-1 + 2 * i / 2000))
                   for i in range(2001)) for d in (1, 3, 6)}
    assert errs[6] < errs[3] < errs[1]


def test_errors():
    with pytest.raises(ValueError):
        floater_hormann_weights([0, 1, 2], d=-1)
    with pytest.raises(ValueError):
        floater_hormann_weights([0, 1, 2], d=20)
