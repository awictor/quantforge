"""Forward-mode automatic differentiation with dual numbers.

A dual number ``a + b*eps`` (with ``eps^2 = 0``) carries a value and its derivative
together: evaluating any function on ``x + 1*eps`` yields ``f(x) + f'(x)*eps`` exactly,
by the chain rule built into the arithmetic. Unlike finite differences there is *no*
truncation or round-off error, and unlike the complex-step trick the function need not
be complex-analytic -- it just has to be composed of the supported operations. Provides
a ``Dual`` type with arithmetic and elementary functions, plus a ``derivative`` helper.
Pure standard library.
"""

import math


class Dual:
    """A dual number ``value + deriv * eps`` for forward-mode autodiff.

    Construct a variable as ``Dual(x, 1.0)`` and a constant as ``Dual(c)``. Arithmetic
    and the module's elementary functions propagate the derivative exactly via the
    chain rule; read ``.value`` and ``.deriv`` from the result.
    """

    __slots__ = ("value", "deriv")

    def __init__(self, value, deriv=0.0):
        self.value = float(value)
        self.deriv = float(deriv)

    @staticmethod
    def _coerce(x):
        return x if isinstance(x, Dual) else Dual(x, 0.0)

    def __add__(self, o):
        o = Dual._coerce(o)
        return Dual(self.value + o.value, self.deriv + o.deriv)
    __radd__ = __add__

    def __sub__(self, o):
        o = Dual._coerce(o)
        return Dual(self.value - o.value, self.deriv - o.deriv)

    def __rsub__(self, o):
        return Dual._coerce(o).__sub__(self)

    def __mul__(self, o):
        o = Dual._coerce(o)
        return Dual(self.value * o.value, self.deriv * o.value + self.value * o.deriv)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = Dual._coerce(o)
        v = o.value
        return Dual(self.value / v,
                    (self.deriv * v - self.value * o.deriv) / (v * v))

    def __rtruediv__(self, o):
        return Dual._coerce(o).__truediv__(self)

    def __neg__(self):
        return Dual(-self.value, -self.deriv)

    def __pow__(self, p):
        if isinstance(p, Dual):
            # d/dx u^v = u^v (v' ln u + v u'/u).
            val = self.value ** p.value
            return Dual(val, val * (p.deriv * math.log(self.value)
                                    + p.value * self.deriv / self.value))
        val = self.value ** p
        return Dual(val, p * self.value ** (p - 1) * self.deriv)

    def __repr__(self):
        return f"Dual({self.value}, {self.deriv})"


def _unary(value_fn, deriv_factor):
    def f(x):
        x = Dual._coerce(x)
        return Dual(value_fn(x.value), deriv_factor(x.value) * x.deriv)
    return f


exp = _unary(math.exp, math.exp)
log = _unary(math.log, lambda v: 1.0 / v)
sqrt = _unary(math.sqrt, lambda v: 0.5 / math.sqrt(v))
sin = _unary(math.sin, math.cos)
cos = _unary(math.cos, lambda v: -math.sin(v))
tan = _unary(math.tan, lambda v: 1.0 / math.cos(v) ** 2)
tanh = _unary(math.tanh, lambda v: 1.0 - math.tanh(v) ** 2)


def derivative(f, x):
    """Exact derivative ``f'(x)`` by forward-mode autodiff.

    ``f`` must accept a :class:`Dual` and return a :class:`Dual`, built from ``Dual``
    arithmetic and this module's elementary functions (``exp``, ``log``, ``sin``, ...).
    Returns ``f'(x)`` with no truncation error.
    """
    return f(Dual(x, 1.0)).deriv
