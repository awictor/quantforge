"""Polynomial and interaction feature expansion.

Expands a feature matrix into all monomials up to a given degree -- squares,
cubes, and cross-products -- so a linear model can fit polynomial and interaction
effects. For ``p`` inputs and degree ``d`` the number of generated terms
(including the bias) is ``C(p + d, d)``; ``interaction_only`` drops pure powers,
keeping only distinct-variable products. Pure standard library.
"""

from itertools import combinations_with_replacement, combinations


def _term_indices(p, degree, interaction_only, include_bias):
    """List of index-tuples, each a monomial (empty tuple = bias)."""
    terms = []
    if include_bias:
        terms.append(())
    for d in range(1, degree + 1):
        if interaction_only:
            # Distinct variables only, so degree cannot exceed p.
            if d <= p:
                terms.extend(combinations(range(p), d))
        else:
            terms.extend(combinations_with_replacement(range(p), d))
    return terms


def polynomial_features(X, degree=2, interaction_only=False, include_bias=True):
    """Expand rows of ``X`` into polynomial / interaction features.

    Parameters
    ----------
    X : list[list[float]]
        Input matrix, ``n`` rows of ``p`` features.
    degree : int
        Maximum total degree of the monomials (>= 1).
    interaction_only : bool
        If True, exclude pure powers (x_i^2, ...), keeping only products of
        distinct variables.
    include_bias : bool
        Prepend a constant 1 column.

    Returns
    -------
    (features, powers) : (list[list[float]], list[tuple])
        The expanded matrix and, for each column, the tuple of input indices whose
        product forms it (``()`` for the bias). With ``degree=1`` and a bias this
        is the original matrix with a leading ones column.
    """
    n = len(X)
    if n == 0:
        raise ValueError("need at least one row")
    p = len(X[0])
    if any(len(r) != p for r in X):
        raise ValueError("all rows must have the same length")
    if degree < 1:
        raise ValueError("degree must be >= 1")

    terms = _term_indices(p, degree, interaction_only, include_bias)
    out = []
    for t in range(n):
        row = []
        for term in terms:
            val = 1.0
            for idx in term:
                val *= X[t][idx]
            row.append(val)
        out.append(row)
    return out, terms
