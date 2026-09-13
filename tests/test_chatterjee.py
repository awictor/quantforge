"""Chatterjee's xi and Blomqvist's beta."""

import math
import random

import pytest

from quantforge import chatterjee_xi, blomqvist_beta, kendall_tau_b


def _xi_ref(x, y):
    n = len(x)
    order = sorted(range(n), key=lambda i: x[i])
    ys = [y[i] for i in order]
    r = [sum(1 for v in ys if v <= yi) for yi in ys]
    l = [sum(1 for v in ys if v >= yi) for yi in ys]
    num = n * sum(abs(r[i + 1] - r[i]) for i in range(n - 1))
    den = 2 * sum(li * (n - li) for li in l)
    return 1 - num / den if den else 0.0


def test_xi_approaches_one_for_monotone():
    # xi -> 1 as n grows for a noiseless monotone relation.
    x = list(range(1000))
    assert chatterjee_xi(x, list(range(1000))) > 0.99


def test_xi_small_under_independence():
    rng = random.Random(3)
    xs = [rng.gauss(0, 1) for _ in range(2000)]
    ys = [rng.gauss(0, 1) for _ in range(2000)]
    assert abs(chatterjee_xi(xs, ys)) < 0.06


def test_xi_detects_nonmonotone_where_kendall_fails():
    rng = random.Random(5)
    x = sorted(rng.uniform(0, 6 * math.pi) for _ in range(2000))
    y = [math.sin(xi) for xi in x]
    assert chatterjee_xi(x, y) > 0.9              # strong functional dependence
    assert abs(kendall_tau_b(x, y)) < 0.3         # rank correlation nearly blind


def test_xi_matches_reference():
    rng = random.Random(9)
    for _ in range(300):
        n = rng.randint(2, 20)
        xs = [rng.gauss(0, 1) for _ in range(n)]
        ys = [rng.gauss(0, 1) for _ in range(n)]
        assert abs(chatterjee_xi(xs, ys) - _xi_ref(xs, ys)) < 1e-9


def test_blomqvist_extremes():
    assert blomqvist_beta([1, 2, 3, 4, 5], [10, 20, 30, 40, 50]) == 1.0
    assert blomqvist_beta([1, 2, 3, 4, 5], [50, 40, 30, 20, 10]) == -1.0


def test_blomqvist_independence():
    rng = random.Random(1)
    xs = [rng.gauss(0, 1) for _ in range(4000)]
    ys = [rng.gauss(0, 1) for _ in range(4000)]
    assert abs(blomqvist_beta(xs, ys)) < 0.05


def test_xi_ties_in_y():
    # A step function: xi should be defined and non-negative.
    x = list(range(20))
    y = [0] * 10 + [1] * 10
    xi = chatterjee_xi(x, y)
    assert 0.0 <= xi <= 1.0


def test_validation():
    with pytest.raises(ValueError):
        chatterjee_xi([1, 2], [1])
    with pytest.raises(ValueError):
        blomqvist_beta([1], [1])
