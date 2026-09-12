"""Levenberg-Marquardt nonlinear least-squares calibration.

Fits parameters ``beta`` so a model ``model(beta, x_i)`` matches targets ``y_i`` in
the least-squares sense. Levenberg-Marquardt interpolates between Gauss-Newton
(fast near the optimum) and gradient descent (robust far from it) via a damping
factor ``lambda`` that grows when a step fails and shrinks when it succeeds:

    (J'J + lambda diag(J'J)) delta = J' r,   beta <- beta - delta.

The Jacobian of the residuals is taken by central finite differences, so only the
model is needed -- no analytic derivatives. The standard calibration engine
(vol-surface, curve, or any parametric fit). Pure standard library.
"""

from .numdiff import jacobian
from .portopt import _invert


def levenberg_marquardt(model, xs, ys, beta0, max_iter=100, tol=1e-10,
                        lam0=1e-3):
    """Fit ``beta`` minimizing ``sum_i (model(beta, xs[i]) - ys[i])^2``.

    Parameters
    ----------
    model : callable
        ``model(beta, x)`` returning a scalar prediction.
    xs, ys : sequences
        Inputs and targets of equal length.
    beta0 : sequence of float
        Initial parameter guess.
    max_iter, tol : int, float
        Iteration cap and convergence tolerance on the parameter step.
    lam0 : float
        Initial damping.

    Returns
    -------
    dict
        ``parameters``, ``residual`` (sum of squared errors), ``iterations``,
        ``converged``.
    """
    n = len(xs)
    if n != len(ys):
        raise ValueError("xs and ys must have the same length")
    p = len(beta0)
    if n < p:
        raise ValueError("need at least as many points as parameters")

    beta = list(map(float, beta0))

    def resid(b):
        return [model(b, xs[i]) - ys[i] for i in range(n)]

    def sse(b):
        r = resid(b)
        return sum(v * v for v in r)

    lam = lam0
    cur_sse = sse(beta)
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        # Residual vector as a function of beta, for the Jacobian.
        J = jacobian(lambda b: resid(b), beta)          # n x p
        r = resid(beta)
        # Normal-equation pieces J'J and J'r.
        JtJ = [[sum(J[t][i] * J[t][j] for t in range(n)) for j in range(p)]
               for i in range(p)]
        Jtr = [sum(J[t][i] * r[t] for t in range(n)) for i in range(p)]

        stepped = False
        for _ in range(30):                              # inner damping search
            A = [[JtJ[i][j] + (lam * JtJ[i][i] if i == j else 0.0)
                  for j in range(p)] for i in range(p)]
            try:
                inv = _invert(A)
            except ValueError:
                lam *= 10.0
                continue
            delta = [sum(inv[i][j] * Jtr[j] for j in range(p)) for i in range(p)]
            trial = [beta[i] - delta[i] for i in range(p)]
            trial_sse = sse(trial)
            if trial_sse < cur_sse:
                step_size = max(abs(d) for d in delta)
                beta, cur_sse = trial, trial_sse
                lam = max(lam / 10.0, 1e-12)
                stepped = True
                if step_size < tol:
                    converged = True
                break
            lam *= 10.0
        if converged or not stepped:
            break

    return {
        "parameters": beta,
        "residual": cur_sse,
        "iterations": it,
        "converged": converged,
    }
