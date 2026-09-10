"""Minimal derivative-free optimizer (Nelder-Mead simplex), pure stdlib.

Used to calibrate model parameters to market quotes without pulling in SciPy.
"""

from typing import Callable, Sequence, List


def nelder_mead(
    f: Callable[[List[float]], float],
    x0: Sequence[float],
    step: float = 0.1,
    max_iter: int = 2000,
    tol: float = 1e-10,
    alpha: float = 1.0,   # reflection
    gamma: float = 2.0,   # expansion
    rho: float = 0.5,     # contraction
    sigma: float = 0.5,   # shrink
):
    """Minimize ``f`` over R^n from ``x0``. Returns (best_x, best_f).

    A textbook Nelder-Mead: build an initial simplex by perturbing each
    coordinate, then reflect/expand/contract/shrink until the spread of
    function values falls below ``tol`` or ``max_iter`` is hit.
    """
    n = len(x0)
    x0 = [float(v) for v in x0]

    # Build initial simplex: x0 plus one vertex per dimension.
    simplex = [x0[:]]
    for i in range(n):
        pt = x0[:]
        delta = step * (abs(pt[i]) if pt[i] != 0 else 1.0)
        pt[i] += delta
        simplex.append(pt)

    fvals = [f(pt) for pt in simplex]

    for _ in range(max_iter):
        # Order vertices by function value.
        order = sorted(range(n + 1), key=lambda i: fvals[i])
        simplex = [simplex[i] for i in order]
        fvals = [fvals[i] for i in order]

        # Convergence: spread of function values.
        if abs(fvals[-1] - fvals[0]) <= tol:
            break

        # Centroid of all but the worst point.
        centroid = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]

        # Reflection.
        worst = simplex[-1]
        reflected = [centroid[j] + alpha * (centroid[j] - worst[j]) for j in range(n)]
        fr = f(reflected)

        if fvals[0] <= fr < fvals[-2]:
            simplex[-1], fvals[-1] = reflected, fr
            continue

        # Expansion.
        if fr < fvals[0]:
            expanded = [centroid[j] + gamma * (reflected[j] - centroid[j]) for j in range(n)]
            fe = f(expanded)
            if fe < fr:
                simplex[-1], fvals[-1] = expanded, fe
            else:
                simplex[-1], fvals[-1] = reflected, fr
            continue

        # Contraction.
        contracted = [centroid[j] + rho * (worst[j] - centroid[j]) for j in range(n)]
        fc = f(contracted)
        if fc < fvals[-1]:
            simplex[-1], fvals[-1] = contracted, fc
            continue

        # Shrink toward the best vertex.
        best = simplex[0]
        for i in range(1, n + 1):
            simplex[i] = [best[j] + sigma * (simplex[i][j] - best[j]) for j in range(n)]
            fvals[i] = f(simplex[i])

    order = sorted(range(n + 1), key=lambda i: fvals[i])
    return simplex[order[0]], fvals[order[0]]
