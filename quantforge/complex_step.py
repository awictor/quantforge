"""Complex-step differentiation: first derivatives with no subtractive cancellation.

A finite difference subtracts two nearly-equal numbers, so shrinking the step to cut
truncation error eventually loses precision to round-off. The complex-step trick
sidesteps that entirely: for a real-analytic ``f`` that accepts a complex argument,

    f(x + i h) = f(x) + i h f'(x) - h^2 f''(x)/2 + ...,

so ``f'(x) = Im(f(x + i h)) / h`` with an error of only ``O(h^2)`` and *no
subtraction* -- take ``h`` as small as ``1e-100`` and the derivative is accurate to
machine precision. The catch: ``f`` must be written with complex-safe operations (no
``abs``, ``max``, or comparisons that break on complex inputs). Pure standard library.
"""


def complex_step_derivative(f, x, h=1e-20):
    """First derivative ``f'(x)`` by the complex-step method.

    ``f`` must accept a complex argument and be analytic near ``x``. Returns the
    derivative ``Im(f(x + i h)) / h``; because there is no subtraction of nearby
    values, the tiny default ``h = 1e-20`` gives essentially machine-precision
    accuracy. Raises ``TypeError`` (propagated) if ``f`` is not complex-safe.
    """
    if h <= 0:
        raise ValueError("step h must be positive")
    return f(complex(x, h)).imag / h


def complex_step_gradient(f, x, h=1e-20):
    """Gradient of a scalar ``f`` of a vector ``x`` by the complex-step method.

    Perturbs each coordinate by ``i h`` in turn. ``f`` must accept a list/sequence of
    (possibly complex) components and return a complex-safe scalar. Returns the list of
    partial derivatives, each to near machine precision.
    """
    n = len(x)
    grad = []
    for j in range(n):
        xc = [complex(x[i], 0.0) for i in range(n)]
        xc[j] = complex(x[j], h)
        grad.append(f(xc).imag / h)
    return grad
