"""Hierarchical Risk Parity and inverse-volatility allocation.

López de Prado's Hierarchical Risk Parity (HRP) builds a diversified portfolio
without inverting the covariance matrix: it clusters assets by correlation
distance, orders them along the cluster tree (quasi-diagonalization), and splits
capital top-down by inverse-variance (recursive bisection). Also includes the
simpler inverse-volatility weighting. Pure standard library.
"""

import math


def inverse_volatility_weights(cov):
    """Inverse-volatility weights ``(1/sigma_i) / sum_j (1/sigma_j)``.

    Each asset weighted by the reciprocal of its standard deviation, normalized to
    sum to one. Higher-volatility assets get less capital; ignores correlations.
    """
    n = len(cov)
    vols = [math.sqrt(cov[i][i]) for i in range(n)]
    if any(v <= 0 for v in vols):
        raise ValueError("variances must be positive")
    inv = [1.0 / v for v in vols]
    s = sum(inv)
    return [x / s for x in inv]


def _corr_from_cov(cov):
    n = len(cov)
    vols = [math.sqrt(cov[i][i]) for i in range(n)]
    return [[cov[i][j] / (vols[i] * vols[j]) for j in range(n)] for i in range(n)]


def _distance(corr):
    """Correlation distance ``d = sqrt((1 - rho)/2)`` in [0, 1]."""
    n = len(corr)
    return [[math.sqrt(max((1.0 - corr[i][j]) / 2.0, 0.0)) for j in range(n)]
            for i in range(n)]


def _quasi_diagonal_order(dist):
    """Single-linkage clustering order (quasi-diagonalization) of the assets.

    Agglomerates the two closest clusters repeatedly (single linkage) and returns
    the leaf order from the resulting tree, so correlated assets sit adjacent.
    """
    n = len(dist)
    clusters = {i: [i] for i in range(n)}
    active = list(range(n))
    next_id = n
    children = {}
    while len(active) > 1:
        best = None
        for a_i in range(len(active)):
            for b_i in range(a_i + 1, len(active)):
                ca, cb = active[a_i], active[b_i]
                # Single-linkage distance between clusters.
                d = min(dist[x][y] for x in clusters[ca] for y in clusters[cb])
                if best is None or d < best[0]:
                    best = (d, ca, cb)
        _, ca, cb = best
        children[next_id] = (ca, cb)
        clusters[next_id] = clusters[ca] + clusters[cb]
        active.remove(ca)
        active.remove(cb)
        active.append(next_id)
        next_id += 1

    def leaves(node):
        if node not in children:
            return [node]
        left, right = children[node]
        return leaves(left) + leaves(right)

    return leaves(active[0])


def _cluster_variance(cov, items):
    """Inverse-variance-weighted variance of a sub-portfolio ``items``."""
    inv = [1.0 / cov[i][i] for i in items]
    s = sum(inv)
    w = [x / s for x in inv]
    var = 0.0
    for a, i in enumerate(items):
        for b, j in enumerate(items):
            var += w[a] * w[b] * cov[i][j]
    return var


def hierarchical_risk_parity(cov):
    """Hierarchical Risk Parity weights (López de Prado).

    Clusters assets by correlation distance, quasi-diagonalizes, and allocates by
    recursive bisection: each split gives the two sub-clusters weights inversely
    proportional to their inverse-variance-weighted variances. Returns positive
    weights summing to one, aligned with the original asset order. Needs no matrix
    inversion, so it is stable for near-singular covariances.
    """
    n = len(cov)
    if n == 0:
        raise ValueError("covariance must be non-empty")
    if any(cov[i][i] <= 0 for i in range(n)):
        raise ValueError("variances must be positive")
    corr = _corr_from_cov(cov)
    order = _quasi_diagonal_order(_distance(corr))
    weights = {i: 1.0 for i in order}
    clusters = [order]
    while clusters:
        new_clusters = []
        for cl in clusters:
            if len(cl) <= 1:
                continue
            half = len(cl) // 2
            left, right = cl[:half], cl[half:]
            v_left = _cluster_variance(cov, left)
            v_right = _cluster_variance(cov, right)
            alpha = 1.0 - v_left / (v_left + v_right)
            for i in left:
                weights[i] *= alpha
            for i in right:
                weights[i] *= (1.0 - alpha)
            new_clusters.extend([left, right])
        clusters = new_clusters
    return [weights[i] for i in range(n)]
