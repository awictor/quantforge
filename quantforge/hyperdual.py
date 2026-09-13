"""Hyperdual numbers: exact first and second derivatives in one evaluation.

A dual number carries a value and one derivative; a hyperdual number
``a + b*e1 + c*e2 + d*e1e2`` (with ``e1^2 = e2^2 = (e1e2)^2 = 0``) carries the value,
*both* first-order slots, and the second derivative in the cross term. Evaluating a
function on a hyperdual seeded with ``b = c = 1`` yields ``f(x)``, ``f'(x)`` (twice) and
``f''(x)`` exactly and simultaneously -- no finite-difference error in either
derivative, unlike the two-step or squared-step schemes that lose precision on the
second derivative. Provides a ``HyperDual`` type and a ``second_derivative`` helper.
Pure standard library.
"""

import math


class HyperDual:
    """Hyperdual number ``f0 + f1 e1 + f2 e2 + f12 e1e2`` for 1st/2nd derivatives.

    Seed a variable as ``HyperDual(x, 1.0, 1.0, 0.0)``; arithmetic and this module's
    elementary functions propagate all four parts exactly. After evaluating ``f``, the
    result's ``f1`` (or ``f2``) is ``f'(x)`` and ``f12`` is ``f''(x)``.
    """

    __slots__ = ("f0", "f1", "f2", "f12")

    def __init__(self, f0, f1=0.0, f2=0.0, f12=0.0):
        self.f0 = float(f0)
        self.f1 = float(f1)
        self.f2 = float(f2)
        self.f12 = float(f12)

    @staticmethod
    def _coerce(x):
        return x if isinstance(x, HyperDual) else HyperDual(x)

    def __add__(self, o):
        o = HyperDual._coerce(o)
        return HyperDual(self.f0 + o.f0, self.f1 + o.f1, self.f2 + o.f2, self.f12 + o.f12)
    __radd__ = __add__

    def __sub__(self, o):
        o = HyperDual._coerce(o)
        return HyperDual(self.f0 - o.f0, self.f1 - o.f1, self.f2 - o.f2, self.f12 - o.f12)

    def __rsub__(self, o):
        return HyperDual._coerce(o).__sub__(self)

    def __mul__(self, o):
        o = HyperDual._coerce(o)
        return HyperDual(
            self.f0 * o.f0,
            self.f1 * o.f0 + self.f0 * o.f1,
            self.f2 * o.f0 + self.f0 * o.f2,
            self.f12 * o.f0 + self.f1 * o.f2 + self.f2 * o.f1 + self.f0 * o.f12,
        )
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = HyperDual._coerce(o)
        return self * o._reciprocal()

    def __rtruediv__(self, o):
        return HyperDual._coerce(o).__mul__(self._reciprocal())

    def _reciprocal(self):
        # 1/f via the chain rule: g(f) with g=1/u, g'=-1/u^2, g''=2/u^3.
        v = self.f0
        return self._chain(1.0 / v, -1.0 / (v * v), 2.0 / (v ** 3))

    def __neg__(self):
        return HyperDual(-self.f0, -self.f1, -self.f2, -self.f12)

    def _chain(self, g, gp, gpp):
        """Apply a scalar function with value ``g``, first ``gp``, second ``gpp``."""
        return HyperDual(
            g,
            gp * self.f1,
            gp * self.f2,
            gp * self.f12 + gpp * self.f1 * self.f2,
        )

    def __pow__(self, p):
        v = self.f0
        return self._chain(v ** p, p * v ** (p - 1), p * (p - 1) * v ** (p - 2))

    def __repr__(self):
        return f"HyperDual({self.f0}, {self.f1}, {self.f2}, {self.f12})"


def _unary(val, d1, d2):
    def f(x):
        x = HyperDual._coerce(x)
        return x._chain(val(x.f0), d1(x.f0), d2(x.f0))
    return f


exp = _unary(math.exp, math.exp, math.exp)
log = _unary(math.log, lambda v: 1.0 / v, lambda v: -1.0 / (v * v))
sqrt = _unary(math.sqrt, lambda v: 0.5 / math.sqrt(v), lambda v: -0.25 / v ** 1.5)
sin = _unary(math.sin, math.cos, lambda v: -math.sin(v))
cos = _unary(math.cos, lambda v: -math.sin(v), lambda v: -math.cos(v))


def derivatives(f, x):
    """Return ``(f(x), f'(x), f''(x))`` exactly by one hyperdual evaluation.

    ``f`` must accept a :class:`HyperDual` built from hyperdual arithmetic and this
    module's elementary functions. No finite-difference error in either derivative.
    """
    r = f(HyperDual(x, 1.0, 1.0, 0.0))
    return r.f0, r.f1, r.f12


def second_derivative(f, x):
    """Exact ``f''(x)`` by hyperdual autodiff."""
    return derivatives(f, x)[2]
