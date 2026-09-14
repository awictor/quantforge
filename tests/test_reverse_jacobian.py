import math

from quantforge import (
    reverse_jacobian,
    reverse_hessian,
    reverse_gradient_vector,
)
from quantforge.numdiff import jacobian as nd_jac, hessian as nd_hess


def close(a, b, tol=1e-5):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_jacobian_vs_numdiff_and_analytic():
    def fvar(v):
        x, y, z = v
        return [x * y, z.sin() + x, (x / y).exp()]

    def fplain(v):
        x, y, z = v
        return [x * y, math.sin(z) + x, math.exp(x / y)]

    pt = [1.3, 2.1, 0.7]
    J = reverse_jacobian(fvar, pt)
    Jn = nd_jac(fplain, pt)
    for i in range(3):
        for j in range(3):
            assert close(J[i][j], Jn[i][j], 1e-4)
    # analytic row 0: grad(x*y) = [y, x, 0]
    assert close(J[0][0], pt[1])
    assert close(J[0][1], pt[0])
    assert close(J[0][2], 0.0)
    # row 1: grad(sin z + x) = [1, 0, cos z]  -- the zero here catches cross-output stale grad
    assert close(J[1][0], 1.0)
    assert close(J[1][1], 0.0)
    assert close(J[1][2], math.cos(pt[2]))


def test_scalar_output_gives_single_row():
    J = reverse_jacobian(lambda v: v[0] * v[1] + v[1], [2.0, 3.0])
    assert len(J) == 1
    assert close(J[0][0], 3.0)
    assert close(J[0][1], 3.0)


def test_hessian_analytic_and_symmetric():
    def g(v):
        x, y, z = v
        return x * x * y + y * y * y + x * z

    pt = [1.5, 2.0, 0.5]
    x, y, z = pt
    H = reverse_hessian(g, pt)
    analytic = [[2 * y, 2 * x, 1], [2 * x, 6 * y, 0], [1, 0, 0]]
    for i in range(3):
        for j in range(3):
            assert close(H[i][j], analytic[i][j], 1e-4)
            assert H[i][j] == H[j][i]


def test_hessian_vs_numdiff():
    def t(v):
        a, b = v
        return (a * b).exp() + a.sin() * b

    def tplain(v):
        a, b = v
        return math.exp(a * b) + math.sin(a) * b

    pt = [0.6, 1.1]
    H = reverse_hessian(t, pt)
    Hn = nd_hess(tplain, pt)
    for i in range(2):
        for j in range(2):
            assert close(H[i][j], Hn[i][j], 1e-3)


def test_gradient_vector_rosenbrock():
    def ros(v):
        return (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2

    a, b = 0.5, 0.7
    gr = reverse_gradient_vector(ros, [a, b])
    assert close(gr[0], -2 * (1 - a) - 400 * a * (b - a * a))
    assert close(gr[1], 200 * (b - a * a))
