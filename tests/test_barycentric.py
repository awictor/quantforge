"""Barycentric Lagrange interpolation."""

import math

import pytest

from quantforge import (barycentric_weights, barycentric_eval,
                        chebyshev_nodes, chebyshev_barycentric_weights, neville)


def test_matches_neville():
    xs = [0, 1, 2, 3, 4]
    ys = [2 * x ** 3 - 3 * x ** 2 + x - 5 for x in xs]
    w = barycentric_weights(xs)
    for xt in (0.5, 1.5, 2.7, 3.9):
        nv, _ = neville(xs, ys, xt)
        assert abs(barycentric_eval(xs, ys, w, xt) - nv) < 1e-9


def test_exact_at_nodes():
    xs = [0, 1, 2, 3]
    ys = [10, 20, 30, 45]
    w = barycentric_weights(xs)
    for i in range(4):
        assert barycentric_eval(xs, ys, w, xs[i]) == ys[i]


def test_chebyshev_closed_form_weights_match():
    n = 8
    cn = chebyshev_nodes(-1, 1, n)
    gw = barycentric_weights(cn)
    cw = chebyshev_barycentric_weights(n)
    ratios = [gw[i] / cw[i] for i in range(n)]
    assert max(ratios) - min(ratios) < 1e-6 * abs(ratios[0])


def test_runge_chebyshev_beats_equispaced():
    runge = lambda x: 1 / (1 + 25 * x * x)
    cn = chebyshev_nodes(-1, 1, 21)
    cw = barycentric_weights(cn)
    cy = [runge(x) for x in cn]
    en = [-1 + 2 * i / 20 for i in range(21)]
    ew = barycentric_weights(en)
    ey = [runge(x) for x in en]
    grid = [i / 100 for i in range(-99, 100)]
    err_c = max(abs(barycentric_eval(cn, cy, cw, x) - runge(x)) for x in grid)
    err_e = max(abs(barycentric_eval(en, ey, ew, x) - runge(x)) for x in grid)
    assert err_c < 0.1
    assert err_e > 1.0                       # equispaced blows up


def test_spectral_convergence_on_exp():
    cn = chebyshev_nodes(0, 2, 25)
    cw = barycentric_weights(cn)
    cy = [math.exp(x) for x in cn]
    grid = [i / 50 for i in range(0, 100)]
    assert max(abs(barycentric_eval(cn, cy, cw, x) - math.exp(x)) for x in grid) < 1e-12


def test_closed_form_weights_same_interpolant():
    cn = chebyshev_nodes(0, 2, 15)
    cy = [math.sin(x) for x in cn]
    v1 = barycentric_eval(cn, cy, barycentric_weights(cn), 0.7)
    v2 = barycentric_eval(cn, cy, chebyshev_barycentric_weights(15), 0.7)
    assert abs(v1 - v2) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        barycentric_weights([1, 1])                  # duplicate
    with pytest.raises(ValueError):
        chebyshev_nodes(0, 1, 1)                      # too few nodes
