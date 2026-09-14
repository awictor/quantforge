import math
import random

from quantforge import Var, reverse_gradient


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_simple_square():
    x = Var(3.0)
    z = x * x
    z.backward()
    assert close(x.grad, 6.0)


def test_analytic_two_var():
    xf, yf = 2.0, 0.5
    g = reverse_gradient(lambda v: (v[0] * v[1]).exp() + v[0].sin(), [xf, yf])
    assert close(g[0], yf * math.exp(xf * yf) + math.cos(xf))
    assert close(g[1], xf * math.exp(xf * yf))


def _F(v):
    a, b, c = v
    return (a * b + (a / c).log() + b.tanh() * c.sqrt() - (a ** 3)).exp() / (1.0 + (b * c).cos())


def _Fval(xs):
    return _F([Var(x) for x in xs]).value


def test_gradient_matches_finite_difference():
    random.seed(3)
    for _ in range(200):
        xs = [random.uniform(0.5, 2.0) for _ in range(3)]
        g = reverse_gradient(_F, xs)
        h = 1e-6
        for i in range(3):
            xp = xs[:]
            xp[i] += h
            xm = xs[:]
            xm[i] -= h
            fd = (_Fval(xp) - _Fval(xm)) / (2 * h)
            assert close(g[i], fd, 1e-4)


def test_shared_subexpression():
    x = Var(2.0)
    s = x.sin()
    z = s * s + s
    z.backward()
    assert close(x.grad, math.cos(2.0) * (2 * math.sin(2.0) + 1))


def test_division_and_reflected_ops():
    g = reverse_gradient(lambda v: (1.0 - v[0]) / (2.0 + v[1]), [0.4, 1.1])
    assert close(g[0], -1 / (2 + 1.1))
    assert close(g[1], -(1 - 0.4) / (2 + 1.1) ** 2)


def test_backward_is_idempotent():
    x = Var(1.5)
    z = x.exp()
    z.backward()
    a = x.grad
    z.backward()
    assert close(x.grad, a)


def test_negative_power():
    x = Var(2.0)
    z = x ** -2
    z.backward()
    assert close(x.grad, -2 * 2.0 ** -3)


def test_elementary_derivatives():
    for xv in (0.3, 1.2, 2.5):
        for meth, dfn in (
            ("exp", math.exp),
            ("sin", lambda t: math.cos(t)),
            ("cos", lambda t: -math.sin(t)),
            ("tanh", lambda t: 1 - math.tanh(t) ** 2),
            ("log", lambda t: 1 / t),
            ("sqrt", lambda t: 0.5 / math.sqrt(t)),
        ):
            x = Var(xv)
            z = getattr(x, meth)()
            z.backward()
            assert close(x.grad, dfn(xv))
