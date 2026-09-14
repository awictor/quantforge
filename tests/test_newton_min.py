import pytest

from quantforge import newton_min, lbfgs


def close(a, b, tol=1e-5):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_quadratic_one_shot():
    target = [3.0, -2.0, 5.0]

    def quad(v):
        return sum((vi - ti) ** 2 for vi, ti in zip(v, target)) + 1.0

    r = newton_min(quad, [0.0, 0.0, 0.0])
    assert r["converged"]
    for xi, ti in zip(r["x"], target):
        assert close(xi, ti)
    assert r["n_iter"] <= 3


def test_rosenbrock():
    def rosen(v):
        return (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2

    r = newton_min(rosen, [-1.2, 1.0], max_iter=200)
    assert close(r["x"][0], 1.0, 1e-4)
    assert close(r["x"][1], 1.0, 1e-4)


def test_agrees_with_lbfgs():
    def f(v):
        return (v[0] - 1) ** 2 + (v[1] + 2) ** 2 + 0.5 * v[0] * v[1] + 0.1 * v[0] ** 4

    rn = newton_min(f, [3.0, 3.0])
    rl = lbfgs(f, [3.0, 3.0])
    for a, b in zip(rn["x"], rl["x"]):
        assert close(a, b, 1e-4)


def test_damping_recovers_from_indefinite_start():
    # f = x^3 - 3x + y^2 : Hessian at x=0 is near-indefinite, local min at (1, 0)
    def h(v):
        return v[0] ** 3 - 3 * v[0] + v[1] ** 2

    r = newton_min(h, [-0.1, 2.0])
    assert r["converged"]
    assert close(r["x"][0], 1.0, 1e-4)
    assert close(r["x"][1], 0.0, 1e-5)


def test_quartic_bowl():
    def q4(v):
        return (v[0] - 2) ** 4 + (v[1] + 1) ** 2

    r = newton_min(q4, [5.0, 5.0])
    assert close(r["x"][0], 2.0, 1e-3)
    assert close(r["x"][1], -1.0, 1e-5)


def test_empty_dimension_raises():
    with pytest.raises(ValueError):
        newton_min(lambda v: 0.0, [])
