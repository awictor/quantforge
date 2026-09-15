"""Bayesian model comparison via the Laplace approximation to the marginal likelihood.

Choosing between models means comparing their *marginal likelihood* (evidence)
``Z = integral p(data | theta) p(theta) dtheta`` -- the probability of the data with the
parameters integrated out, which automatically penalizes complexity (Occam's razor). The
integral is rarely closed-form, but the Laplace approximation is: expand the log-posterior to
second order about its mode ``theta*`` (a Gaussian), giving

    log Z ~ log p(data | theta*) + log p(theta*) + (d/2) log(2 pi) - (1/2) log det(H)

where ``H`` is the Hessian of the *negative* log-posterior at the mode. The mode is found by
:func:`quantforge.newton_min` and ``H`` by reverse-mode autodiff, so the log-joint is written
once with :class:`quantforge.reverse_ad.Var`. From two evidences come the Bayes factor and
posterior model probabilities. Pure standard library.
"""

import math

from .newton_min import newton_min
from .reverse_jacobian import reverse_hessian
from .reverse_ad import Var
from .lu import lu_solve


def _logdet_via_lu(H):
    # log|det H| from the LU factorization implied by solving against the identity is awkward;
    # instead do a plain Gaussian elimination tracking the product of pivots.
    n = len(H)
    A = [row[:] for row in H]
    logdet = 0.0
    sign = 1.0
    for k in range(n):
        # partial pivot
        piv = max(range(k, n), key=lambda i: abs(A[i][k]))
        if abs(A[piv][k]) < 1e-300:
            return float("-inf")
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        logdet += math.log(abs(A[k][k]))
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] -= f * A[k][j]
    return logdet


def laplace_log_evidence(log_joint, theta0, h=1e-5):
    """Laplace approximation to ``log Z = log integral exp(log_joint(theta)) dtheta``.

    ``log_joint`` maps a length-``d`` list of :class:`Var` to the (unnormalized) log-joint
    ``log p(data | theta) + log p(theta)`` as a single ``Var``. Finds the posterior mode by
    Newton minimization of ``-log_joint``, evaluates the Hessian of ``-log_joint`` there by
    autodiff, and returns ``log_joint(mode) + (d/2) log(2 pi) - (1/2) log det H`` together with
    the ``mode``.
    """
    d = len(theta0)
    neg = lambda th: -log_joint(th)
    res = newton_min(neg, theta0, h=h)
    mode = res["x"]
    # value of the log-joint at the mode
    peak = log_joint([Var(m) for m in mode]).value
    # Hessian of the negative log-joint at the mode (positive-definite at a max of log_joint)
    H = reverse_hessian(neg, mode, h)
    logdet = _logdet_via_lu(H)
    log_z = peak + 0.5 * d * math.log(2.0 * math.pi) - 0.5 * logdet
    return {"log_evidence": log_z, "mode": mode, "log_joint_at_mode": peak}


def bayes_factor(log_evidence_1, log_evidence_2):
    """Bayes factor ``B_12 = Z_1 / Z_2 = exp(log_evidence_1 - log_evidence_2)``.

    ``> 1`` favours model 1, ``< 1`` favours model 2. Computed from log-evidences to avoid
    overflow.
    """
    return math.exp(log_evidence_1 - log_evidence_2)


def posterior_model_probabilities(log_evidences, priors=None):
    """Posterior model probabilities from a list of log-evidences (equal priors by default).

    Combines ``log p(M_k) + log Z_k`` and normalizes stably via the log-sum-exp trick. Returns a
    list of probabilities summing to 1.
    """
    m = len(log_evidences)
    if priors is None:
        log_prior = [0.0] * m
    else:
        s = sum(priors)
        log_prior = [math.log(p / s) if p > 0 else float("-inf") for p in priors]
    logpost = [log_evidences[k] + log_prior[k] for k in range(m)]
    mx = max(logpost)
    unnorm = [math.exp(lp - mx) for lp in logpost]
    total = sum(unnorm)
    return [u / total for u in unnorm]
