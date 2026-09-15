"""Sliced Wasserstein distance for multivariate samples.

The Sinkhorn solver (:mod:`quantforge.sinkhorn`) gives entropic optimal transport for arbitrary
cost matrices, but its cost grows with the number of points squared. The *sliced* Wasserstein
distance is a cheaper alternative for point clouds: it projects both samples onto many random
1-D directions, computes the exact 1-D Wasserstein distance along each (a sort), and averages.
This inherits the metric properties of the true Wasserstein distance while costing only
``O(L n log n)`` for ``L`` projections. It is the standard distance for generative-model
training and shape comparison.

A seeded :class:`quantforge.pcg.PCG32` stream makes the random projections reproducible; the
Gaussian Box-Muller helper is shared with the particle filter. Pure standard library.
"""

import math

from .pcg import PCG32
from .particle_filter import pcg_gaussian
from .wasserstein import wasserstein_distance


def _random_unit_vector(dim, rng):
    v = [pcg_gaussian(rng, 0.0, 1.0) for _ in range(dim)]
    norm = math.sqrt(sum(c * c for c in v))
    if norm == 0.0:
        v[0] = 1.0
        norm = 1.0
    return [c / norm for c in v]


def sliced_wasserstein(xs, ys, n_projections=200, p=2, seed=12345):
    """Sliced ``p``-Wasserstein distance between two multivariate point clouds.

    ``xs`` and ``ys`` are lists of points (each a length-``d`` coordinate list; ``d = 1`` scalars
    are also accepted as bare floats). Projects both onto ``n_projections`` random unit
    directions, takes the 1-D ``p``-Wasserstein distance along each, and returns the ``p``-mean
    ``(mean_l W_p(proj_l)^p)^{1/p}``. A seeded PCG32 stream makes the projections reproducible.
    """
    if not xs or not ys:
        raise ValueError("both samples must be non-empty")
    # normalize scalars to 1-D vectors
    if not isinstance(xs[0], (list, tuple)):
        xs = [[v] for v in xs]
        ys = [[v] for v in ys]
    dim = len(xs[0])
    rng = PCG32(seed)
    acc = 0.0
    for _ in range(n_projections):
        u = _random_unit_vector(dim, rng)
        px = [sum(pt[k] * u[k] for k in range(dim)) for pt in xs]
        py = [sum(pt[k] * u[k] for k in range(dim)) for pt in ys]
        acc += wasserstein_distance(px, py, p=p) ** p
    return (acc / n_projections) ** (1.0 / p)
