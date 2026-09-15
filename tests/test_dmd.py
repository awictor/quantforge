import cmath
import math
import random

import pytest

from quantforge import dmd
from quantforge.lu import lu_solve


def close(a, b, tol=1e-4):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _evolve(A, x0, steps):
    x = list(x0)
    snaps = [x[:]]
    d = len(A)
    for _ in range(steps):
        x = [sum(A[i][k] * x[k] for k in range(d)) for i in range(d)]
        snaps.append(x[:])
    return snaps


def test_linear_system_eigenvalues():
    A = [[0.9, 0.1], [0.0, 0.8]]
    res = dmd(_evolve(A, [1.0, 1.0], 15))
    ev = sorted(z.real for z in res["eigenvalues"])
    assert close(ev[0], 0.8) and close(ev[1], 0.9)


def test_rotation_with_decay():
    w, r = 0.3, 0.95
    A = [[r * math.cos(w), -r * math.sin(w)], [r * math.sin(w), r * math.cos(w)]]
    res = dmd(_evolve(A, [1.0, 0.0], 30))
    for z in res["eigenvalues"]:
        assert close(abs(complex(z)), r, 1e-3)
        assert close(abs(cmath.phase(complex(z))), w, 1e-3)
    assert all(close(g, r, 1e-3) for g in res["growth_rates"])


def test_three_mode_system():
    random.seed(0)
    lam = [0.99, 0.7, 0.4]
    P = [[random.gauss(0, 1) for _ in range(3)] for _ in range(3)]

    def inv(M):
        n = len(M)
        cols = [lu_solve(M, [1.0 if i == k else 0.0 for i in range(n)]) for k in range(n)]
        return [[cols[j][i] for j in range(n)] for i in range(n)]

    def mm(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(len(Y))) for j in range(len(Y[0]))]
                for i in range(len(X))]

    D = [[lam[i] if i == j else 0.0 for j in range(3)] for i in range(3)]
    A = mm(mm(P, D), inv(P))
    res = dmd(_evolve(A, [1.0, 0.5, -0.3], 20))
    ev = sorted(z.real for z in res["eigenvalues"])
    for got, want in zip(ev, sorted(lam)):
        assert close(got, want, 1e-3)


def test_too_few_snapshots_raises():
    with pytest.raises(ValueError):
        dmd([[1.0, 2.0]])
