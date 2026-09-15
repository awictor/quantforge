"""Bayesian optimization: minimize an expensive black box with a GP surrogate.

When each evaluation of ``f`` is costly (a simulation, a backtest, a physical experiment), you
cannot afford a dense grid search. Bayesian optimization fits a Gaussian process to the points
seen so far (:mod:`quantforge.gaussian_process`) and uses an *acquisition function* to pick the
next point -- balancing exploitation (where the surrogate predicts a low value) against
exploration (where it is most uncertain). Expected improvement is the classic choice.

* :func:`expected_improvement` -- EI of a candidate given the GP posterior and the best value so
  far, in closed form via the normal CDF/PDF.
* :func:`bayesian_optimize` -- the loop: fit GP, maximize EI over a candidate set, evaluate ``f``
  there, repeat. Minimizes by default.

Pure standard library on top of the GP module and :func:`quantforge.mathfns.norm_cdf`.
"""

import math

from .gaussian_process import gp_predict
from .mathfns import norm_cdf


def _norm_pdf(z):
    return math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi)


def expected_improvement(mean, var, best, xi=0.01):
    """Expected improvement of a candidate for a *minimization* problem.

    ``mean``/``var`` are the GP posterior at the candidate, ``best`` the lowest observed value,
    ``xi`` an exploration margin. ``EI = (best - mean - xi) Phi(z) + sigma phi(z)`` with
    ``z = (best - mean - xi) / sigma``; zero when the posterior is certain (``var = 0``).
    """
    sigma = math.sqrt(max(var, 0.0))
    if sigma < 1e-12:
        return 0.0
    imp = best - mean - xi
    z = imp / sigma
    return imp * norm_cdf(z) + sigma * _norm_pdf(z)


def bayesian_optimize(f, candidates, n_init=3, n_iter=15, length_scale=1.0,
                      variance=1.0, noise=1e-6, xi=0.01, seed=12345):
    """Minimize ``f`` over a discrete ``candidates`` set by GP-based Bayesian optimization.

    Evaluates ``n_init`` seed points (evenly spread across the candidate list, deterministic),
    then for ``n_iter`` rounds fits a GP to all evaluations and picks the unevaluated candidate
    with the highest :func:`expected_improvement`. ``candidates`` is a list of inputs (scalars or
    coordinate lists) accepted by ``f`` and the GP. Returns a dict with ``best_x``, ``best_y``,
    the full ``X``/``y`` history, and ``n_eval``.
    """
    m = len(candidates)
    if m == 0:
        raise ValueError("need at least one candidate")
    # deterministic spread of initial indices
    init_idx = sorted(set(int(round(i * (m - 1) / max(1, n_init - 1))) for i in range(n_init)))
    evaluated = {}
    for i in init_idx:
        evaluated[i] = f(candidates[i])

    for _ in range(n_iter):
        X = [candidates[i] for i in evaluated]
        y = [evaluated[i] for i in evaluated]
        best = min(y)
        # remaining candidates
        remaining = [i for i in range(m) if i not in evaluated]
        if not remaining:
            break
        pred = gp_predict(X, y, [candidates[i] for i in remaining],
                          length_scale=length_scale, variance=variance, noise=noise)
        # pick max EI
        best_ei = -1.0
        pick = remaining[0]
        for k, i in enumerate(remaining):
            ei = expected_improvement(pred["mean"][k], pred["var"][k], best, xi)
            if ei > best_ei:
                best_ei = ei
                pick = i
        evaluated[pick] = f(candidates[pick])

    best_i = min(evaluated, key=lambda i: evaluated[i])
    return {
        "best_x": candidates[best_i],
        "best_y": evaluated[best_i],
        "X": [candidates[i] for i in evaluated],
        "y": [evaluated[i] for i in evaluated],
        "n_eval": len(evaluated),
    }
