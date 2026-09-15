"""Observability and Luenberger observer design (the dual of controllability/control).

A state-feedback controller needs the full state ``x``; a real system only measures outputs
``y = C x``. An *observer* reconstructs ``x`` from ``y`` -- it is the deterministic cousin of the
Kalman filter. Whether that is possible is governed by *observability*, the exact dual of
controllability: the pair ``(A, C)`` is observable iff ``(A^T, C^T)`` is controllable. By that
duality the observer gain ``L`` for ``x_hat_{k+1} = A x_hat + L(y - C x_hat)`` is a pole-placement
problem on the transposed system, so the error dynamics ``A - L C`` get the eigenvalues you
choose.

:func:`observability_matrix` / :func:`is_observable` mirror the controllability tools;
:func:`observer_gain` places the estimator poles via Ackermann on the dual. Pure standard library.
"""

from .lqr import _matmul, _transpose  # reuse matrix helpers
from .pole_placement import ackermann


def observability_matrix(A, C):
    """Observability matrix ``[C; CA; CA^2; ...; CA^{n-1}]`` (stacked vertically, ``n*p x n``)."""
    n = len(A)
    rows = []
    CAk = [row[:] for row in C]
    rows.extend([r[:] for r in CAk])
    for _ in range(1, n):
        CAk = _matmul(CAk, A)
        rows.extend([r[:] for r in CAk])
    return rows


def is_observable(A, C, tol=1e-9):
    """True if ``(A, C)`` is observable (observability matrix has full column rank ``n``)."""
    from .svd import matrix_rank
    return matrix_rank(observability_matrix(A, C), tol) == len(A)


def observer_gain(A, C, desired_poles):
    """Luenberger observer gain ``L`` placing the estimator error poles at ``desired_poles``.

    Single-output systems (``C`` is ``1 x n``). By duality ``(A - L C)`` and ``(A^T - C^T K^T)``
    share eigenvalues, so ``L = K^T`` where ``K`` places the poles of the transposed system.
    Returns ``L`` as an ``n x 1`` column. Requires ``(A, C)`` observable.
    """
    if len(C) != 1:
        raise ValueError("observer_gain is for single-output systems (C is 1 x n)")
    At = _transpose(A)
    Ct = _transpose(C)               # n x 1, the "B" of the dual system
    K = ackermann(At, Ct, desired_poles)   # 1 x n
    return _transpose(K)             # n x 1
