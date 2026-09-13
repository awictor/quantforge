"""Portfolio concentration and diversification measures.

Three lenses on how spread out a portfolio really is:

- ``herfindahl_index`` -- ``sum w_i^2`` of the weights; 1/n at equal weight, 1 when
  fully concentrated,
- ``effective_number_of_constituents`` -- ``1 / HHI``, the equivalent count of
  equally-weighted holdings,
- ``effective_number_of_bets`` -- Meucci's diversification measure: the exponential
  of the entropy of the risk contributions of the *uncorrelated* principal-component
  factors, which counts truly independent risk sources rather than nominal
  positions.

Pure standard library on top of the Jacobi eigensolver.
"""

import math

from .pca import jacobi_eigen


def herfindahl_index(weights):
    """Herfindahl-Hirschman concentration index ``sum w_i^2``.

    Weights need not be normalized; they are normalized to sum to one first
    (absolute values, for long-short books). Ranges from ``1/n`` (equal weight) to
    ``1`` (a single holding).
    """
    if not weights:
        raise ValueError("weights must be non-empty")
    total = sum(abs(w) for w in weights)
    if total <= 0.0:
        raise ValueError("weights must not all be zero")
    p = [abs(w) / total for w in weights]
    return sum(pi * pi for pi in p)


def effective_number_of_constituents(weights):
    """Effective number of holdings ``1 / HHI``.

    The count of equally-weighted positions with the same concentration. Equals the
    number of holdings at equal weight and approaches 1 as weight concentrates.
    """
    return 1.0 / herfindahl_index(weights)


def effective_number_of_bets(weights, cov):
    """Meucci's effective number of bets from the weights and covariance.

    Diagonalizes the covariance into uncorrelated principal-component factors, splits
    the portfolio variance into each factor's contribution ``p_k`` (summing to one),
    and returns ``exp(-sum p_k ln p_k)`` -- the exponential of the entropy of those
    contributions. Counts independent risk sources: ``n`` when the variance is spread
    evenly across uncorrelated factors, and down toward 1 when one factor dominates.
    """
    n = len(weights)
    if n == 0 or len(cov) != n or any(len(row) != n for row in cov):
        raise ValueError("cov must be an n x n matrix matching the weights")
    eigenvalues, eigenvectors = jacobi_eigen(cov)
    # eigenvectors[i] is the unit eigenvector for eigenvalues[i]. Exposure of the
    # portfolio to factor i is w . v_i; its variance contribution is lambda_i (w.v_i)^2.
    contribs = []
    for i in range(n):
        v = eigenvectors[i]
        exposure = sum(weights[j] * v[j] for j in range(n))
        contribs.append(max(eigenvalues[i], 0.0) * exposure * exposure)
    total = sum(contribs)
    if total <= 0.0:
        raise ValueError("portfolio variance must be positive")
    p = [c / total for c in contribs]
    entropy = -sum(pi * math.log(pi) for pi in p if pi > 0.0)
    return math.exp(entropy)
