"""Discrete state-space (LTI system) simulation and responses.

A discrete linear time-invariant system ``x_{k+1} = A x_k + B u_k``, ``y_k = C x_k + D u_k`` is
the model behind the controllers (:mod:`quantforge.lqr`) and estimators (:mod:`quantforge.observability`).
This simulates it: :func:`lti_simulate` propagates an arbitrary input sequence,
:func:`lti_step_response` and :func:`lti_impulse_response` give the standard test responses, and
:func:`dc_gain` is the steady-state output for a unit step. These characterize a system's
behaviour -- settling time, overshoot, static gain -- from its ``(A, B, C, D)`` matrices.

Single-input single-output or multivariable (matrices sized accordingly). Pure standard library.
"""

from .lu import lu_solve


def _matvec(M, v):
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


def lti_simulate(A, B, C, D, u, x0=None):
    """Simulate ``x_{k+1}=Ax+Bu, y_k=Cx+Du`` over an input sequence ``u``.

    ``A`` (n x n), ``B`` (n x m), ``C`` (p x n), ``D`` (p x m). ``u`` is a list of length-``m``
    input vectors (or scalars for a single input). Returns the list of output vectors ``y_k`` (one
    per input sample). ``x0`` defaults to zero.
    """
    n = len(A)
    x = [0.0] * n if x0 is None else [float(v) for v in x0]
    ys = []
    for uk in u:
        uv = uk if isinstance(uk, (list, tuple)) else [uk]
        y = [sum(C[i][j] * x[j] for j in range(n)) + sum(D[i][j] * uv[j] for j in range(len(uv)))
             for i in range(len(C))]
        ys.append(y)
        Ax = _matvec(A, x)
        Bu = [sum(B[i][j] * uv[j] for j in range(len(uv))) for i in range(n)]
        x = [Ax[i] + Bu[i] for i in range(n)]
    return ys


def lti_step_response(A, B, C, D, n_steps, input_channel=0):
    """Step response: outputs to a unit step on ``input_channel`` for ``n_steps`` samples."""
    m = len(B[0])
    u = [[1.0 if j == input_channel else 0.0 for j in range(m)] for _ in range(n_steps)]
    return lti_simulate(A, B, C, D, u)


def lti_impulse_response(A, B, C, D, n_steps, input_channel=0):
    """Impulse response: outputs to a unit impulse on ``input_channel`` for ``n_steps`` samples."""
    m = len(B[0])
    u = [[0.0] * m for _ in range(n_steps)]
    u[0][input_channel] = 1.0
    return lti_simulate(A, B, C, D, u)


def dc_gain(A, B, C, D):
    """Steady-state (DC) gain matrix ``C (I - A)^{-1} B + D`` -- the output for a unit step.

    Requires ``A`` stable (``I - A`` invertible). Returns the ``p x m`` gain.
    """
    n = len(A)
    m = len(B[0])
    ImA = [[(1.0 if i == j else 0.0) - A[i][j] for j in range(n)] for i in range(n)]
    # (I - A)^{-1} B, column by column
    W = [[0.0] * m for _ in range(n)]
    for j in range(m):
        col = lu_solve(ImA, [B[i][j] for i in range(n)])
        for i in range(n):
            W[i][j] = col[i]
    p = len(C)
    return [[sum(C[i][k] * W[k][j] for k in range(n)) + D[i][j] for j in range(m)]
            for i in range(p)]
