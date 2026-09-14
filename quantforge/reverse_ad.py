"""Reverse-mode automatic differentiation (a small backpropagation tape).

Forward-mode autodiff (:mod:`quantforge.dual`, :mod:`quantforge.hyperdual`) propagates one
derivative per pass, which is cheap when there is one input and many outputs. Reverse mode is
the opposite: one backward pass yields the gradient with respect to *all* inputs at once, which
is why it is the engine behind neural-network training. This module records each elementary
operation on a tape of :class:`Var` nodes, then walks the tape in reverse (topological order)
accumulating adjoints via the chain rule.

Usage::

    x = Var(3.0); y = Var(4.0)
    z = (x * y + x.sin()).exp()
    z.backward()
    x.grad, y.grad          # d z / d x, d z / d y

Pure standard library.
"""

import math


class Var:
    """A scalar node on the autodiff tape.

    Wrap each independent input in a ``Var``; build an expression with the arithmetic operators
    and the elementary methods (``exp``, ``log``, ``sin``, ``cos``, ``tanh``, ``sqrt``, ...).
    Calling :meth:`backward` on the output populates ``.grad`` on every node reachable from it.
    """

    __slots__ = ("value", "grad", "_parents", "_backward")

    def __init__(self, value, _parents=(), _backward=None):
        self.value = float(value)
        self.grad = 0.0
        self._parents = _parents
        self._backward = _backward  # closure that pushes this node's grad to its parents

    def __repr__(self):
        return f"Var(value={self.value}, grad={self.grad})"

    # ---- helpers -------------------------------------------------------
    @staticmethod
    def _coerce(o):
        return o if isinstance(o, Var) else Var(o)

    # ---- arithmetic ----------------------------------------------------
    def __add__(self, o):
        o = self._coerce(o)
        out = Var(self.value + o.value, (self, o))

        def back(g):
            self.grad += g
            o.grad += g

        out._backward = back
        return out

    __radd__ = __add__

    def __mul__(self, o):
        o = self._coerce(o)
        out = Var(self.value * o.value, (self, o))

        def back(g):
            self.grad += o.value * g
            o.grad += self.value * g

        out._backward = back
        return out

    __rmul__ = __mul__

    def __neg__(self):
        out = Var(-self.value, (self,))

        def back(g):
            self.grad += -g

        out._backward = back
        return out

    def __sub__(self, o):
        return self + (-self._coerce(o))

    def __rsub__(self, o):
        return self._coerce(o) + (-self)

    def __truediv__(self, o):
        o = self._coerce(o)
        out = Var(self.value / o.value, (self, o))

        def back(g):
            self.grad += g / o.value
            o.grad += -self.value / (o.value * o.value) * g

        out._backward = back
        return out

    def __rtruediv__(self, o):
        return self._coerce(o) / self

    def __pow__(self, p):
        # p is a constant exponent (float/int)
        out = Var(self.value ** p, (self,))

        def back(g):
            self.grad += p * self.value ** (p - 1) * g

        out._backward = back
        return out

    # ---- elementary functions -----------------------------------------
    def exp(self):
        e = math.exp(self.value)
        out = Var(e, (self,))

        def back(g):
            self.grad += e * g

        out._backward = back
        return out

    def log(self):
        out = Var(math.log(self.value), (self,))

        def back(g):
            self.grad += g / self.value

        out._backward = back
        return out

    def sin(self):
        out = Var(math.sin(self.value), (self,))
        c = math.cos(self.value)

        def back(g):
            self.grad += c * g

        out._backward = back
        return out

    def cos(self):
        out = Var(math.cos(self.value), (self,))
        s = math.sin(self.value)

        def back(g):
            self.grad += -s * g

        out._backward = back
        return out

    def tanh(self):
        t = math.tanh(self.value)
        out = Var(t, (self,))

        def back(g):
            self.grad += (1.0 - t * t) * g

        out._backward = back
        return out

    def sqrt(self):
        r = math.sqrt(self.value)
        out = Var(r, (self,))

        def back(g):
            self.grad += g / (2.0 * r)

        out._backward = back
        return out

    # ---- reverse pass --------------------------------------------------
    def backward(self):
        """Accumulate gradients into every ancestor's ``.grad`` (seed ``d self/d self = 1``)."""
        topo = []
        seen = set()

        def build(node):
            if id(node) in seen:
                return
            seen.add(id(node))
            for p in node._parents:
                build(p)
            topo.append(node)

        build(self)
        # reset grads on the reachable subgraph, then seed and propagate
        for node in topo:
            node.grad = 0.0
        self.grad = 1.0
        for node in reversed(topo):
            if node._backward is not None:
                node._backward(node.grad)


def reverse_gradient(f, xs):
    """Gradient of scalar ``f`` at the point ``xs`` (a sequence of floats).

    ``f`` takes a list of :class:`Var` and returns a single :class:`Var`. Returns the list of
    partial derivatives ``[df/dx_0, ...]`` via one reverse pass.
    """
    vars_ = [Var(x) for x in xs]
    out = f(vars_)
    out.backward()
    return [v.grad for v in vars_]
