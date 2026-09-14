"""Kabsch algorithm: the optimal rotation aligning two point sets.

Given two sets of paired points ``P`` and ``Q``, the Kabsch algorithm finds the rotation
that best superimposes ``P`` onto ``Q`` in the least-squares sense (minimum RMSD) -- the
core of structure comparison, point-cloud registration, and rigid-body fitting. Both sets
are centred, their cross-covariance is SVD-decomposed, and the rotation follows, with a
reflection correction so the result is a proper rotation (determinant +1). Returns the
rotation matrix, the translation, and the residual RMSD. Pure standard library.
"""

from .svd import svd


def _centroid(pts):
    n = len(pts)
    dim = len(pts[0])
    return [sum(p[k] for p in pts) / n for k in range(dim)]


def _det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def kabsch(P, Q):
    """Optimal rotation/translation mapping points ``P`` onto ``Q`` (least-squares).

    ``P`` and ``Q`` are equal-length lists of 3-D points (paired). Returns
    ``(R, t, rmsd)``: the 3x3 rotation matrix and translation ``t`` such that
    ``R @ (p - centroid_P) + centroid_Q`` best matches ``q``, and the residual RMSD.
    """
    if len(P) != len(Q):
        raise ValueError("P and Q must have the same number of points")
    if len(P) == 0:
        raise ValueError("need at least one point")
    if any(len(p) != 3 for p in P) or any(len(q) != 3 for q in Q):
        raise ValueError("points must be 3-D")
    n = len(P)
    cP = _centroid(P)
    cQ = _centroid(Q)
    Pc = [[p[k] - cP[k] for k in range(3)] for p in P]
    Qc = [[q[k] - cQ[k] for k in range(3)] for q in Q]
    # cross-covariance H = Pc^T Qc (3x3)
    H = [[sum(Pc[i][a] * Qc[i][b] for i in range(n)) for b in range(3)] for a in range(3)]
    U, s, V = svd(H)                       # H = U diag(s) V'
    # rotation R = V U'
    R = [[sum(V[i][k] * U[j][k] for k in range(3)) for j in range(3)] for i in range(3)]
    if _det3(R) < 0:
        # reflection: flip the sign of the column of V tied to the smallest singular value
        V2 = [row[:] for row in V]
        for i in range(3):
            V2[i][2] = -V2[i][2]
        R = [[sum(V2[i][k] * U[j][k] for k in range(3)) for j in range(3)] for i in range(3)]
    # translation and residual
    t = [cQ[k] - sum(R[k][j] * cP[j] for j in range(3)) for k in range(3)]
    sq = 0.0
    for i in range(n):
        pred = [sum(R[k][j] * Pc[i][j] for j in range(3)) + cQ[k] for k in range(3)]
        sq += sum((pred[k] - Q[i][k]) ** 2 for k in range(3))
    rmsd = (sq / n) ** 0.5
    return R, t, rmsd
