"""Gauss-Laguerre quadrature for integrals over the half-line.

Approximates integrals against the exponential weight on ``[0, inf)``,

    integral_0^inf f(x) e^{-x} dx ~ sum_i w_i f(x_i),

with ``n`` nodes exact when ``f`` is a polynomial of degree up to ``2n - 1``. A
convenience wrapper handles a general integrand ``integral_0^inf g(x) dx`` by
factoring out ``e^{-x}`` (evaluating ``g(x) e^{x}`` at the nodes) and an
exponential rate. The nodes and weights come from the Golub-Welsch algorithm: the
nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix of the
Laguerre recurrence (diagonal ``2k + 1``, off-diagonal ``k``), and the weights are
the squared first components of the normalized eigenvectors. Pure standard library
on top of :mod:`quantforge.pca`.
"""

import math

from .pca import jacobi_eigen


def gauss_laguerre_nodes_weights(n):
    """Gauss-Laguerre nodes and weights for the weight ``e^{-x}`` on ``[0, inf)``.

    Returns ``(nodes, weights)`` with ``sum_i w_i = 1`` and
    ``sum_i w_i x_i^m = m!`` (the moments of the ``e^{-x}`` density). Exact for
    polynomials up to degree ``2n - 1``. Nodes are positive and returned in
    increasing order.
    """
    if n < 1:
        raise ValueError("n must be at least 1")
    if n == 1:
        return [1.0], [1.0]
    jac = [[0.0] * n for _ in range(n)]
    for k in range(n):
        jac[k][k] = 2.0 * k + 1.0
    for k in range(1, n):
        jac[k][k - 1] = float(k)
        jac[k - 1][k] = float(k)
    eigenvalues, eigenvectors = jacobi_eigen(jac)
    pairs = []
    for i in range(n):
        vec = eigenvectors[i]
        norm = math.sqrt(sum(c * c for c in vec))
        first = vec[0] / norm
        pairs.append((eigenvalues[i], first * first))
    pairs.sort(key=lambda p: p[0])
    nodes = [p[0] for p in pairs]
    weights = [p[1] for p in pairs]
    # Renormalize so the weights sum to exactly one (mu_0 = 1).
    total = sum(weights)
    weights = [w / total for w in weights]
    return nodes, weights


def gauss_laguerre_integral(g, n=32, rate=1.0):
    """Approximate ``integral_0^inf g(x) dx`` by Gauss-Laguerre quadrature.

    Writes the integrand as ``g(x) = [g(x) e^{rate * x}] e^{-rate * x}`` and applies
    the ``e^{-x}`` rule after the substitution ``u = rate * x``:

        integral_0^inf g(x) dx = (1/rate) sum_i w_i g(x_i / rate) e^{x_i}.

    ``rate`` should roughly match the integrand's exponential decay for best
    accuracy (choose ``rate`` near the true decay constant). Exact when
    ``g(x) e^{rate x}`` is a polynomial of degree up to ``2n - 1``.
    """
    if rate <= 0.0:
        raise ValueError("rate must be positive")
    nodes, weights = gauss_laguerre_nodes_weights(n)
    total = 0.0
    for i in range(len(nodes)):
        x = nodes[i] / rate
        total += weights[i] * g(x) * math.exp(nodes[i])
    return total / rate
