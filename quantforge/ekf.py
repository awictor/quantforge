"""Extended Kalman filter (EKF) with exact autodiff Jacobians.

The linear :func:`quantforge.kalman_filter.kalman_filter` assumes ``x_t = F x_{t-1}`` and
``y_t = H x_t``. Many real systems are nonlinear: ``x_t = f(x_{t-1})``, ``y_t = h(x_t)`` for
smooth ``f, h`` (orbital mechanics, tracking with range/bearing sensors, chemical kinetics).
The EKF linearizes ``f`` and ``h`` around the current estimate at each step, using the
Jacobians ``F = df/dx`` and ``H = dh/dx``, then applies the ordinary Kalman predict/update.

Here the Jacobians are *exact*, obtained by reverse-mode autodiff -- ``f`` and ``h`` are written
once with :class:`quantforge.reverse_ad.Var` arithmetic, so there is no hand-coded linearization
and no finite-difference error. Reduces exactly to the linear Kalman filter when ``f, h`` are
linear. Pure standard library.
"""

from .reverse_jacobian import reverse_jacobian
from .reverse_ad import Var


def _matvec(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def _add(A, B, sign=1.0):
    return [[A[i][j] + sign * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _eye(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def _inv(A):
    from .lu import lu_solve
    n = len(A)
    cols = [lu_solve(A, [1.0 if i == j else 0.0 for i in range(n)]) for j in range(n)]
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def extended_kalman_filter(observations, f, h, Q, R, x0, P0):
    """Extended Kalman filter over ``observations`` for nonlinear ``f`` and ``h``.

    ``f`` maps a length-``n`` list of :class:`Var` (the state) to a length-``n`` list -- the
    deterministic state transition. ``h`` maps the state to a length-``m`` list -- the
    measurement model. ``Q`` (n x n), ``R`` (m x m) are the process and observation covariances,
    ``x0`` (n) and ``P0`` (n x n) the initial mean and covariance. Returns a dict with
    ``filtered_means`` and ``filtered_covariances``.

    The transition/observation Jacobians are computed exactly at each step by reverse-mode
    autodiff, so no analytic linearization is required.
    """
    n = len(x0)
    x = list(x0)
    P = [row[:] for row in P0]
    I = _eye(n)

    means, covs = [], []
    for z in observations:
        # ---- predict ----
        Fj = reverse_jacobian(f, x)                     # df/dx at current estimate
        xp = [v.value for v in f([Var(xi) for xi in x])]
        Pp = _add(_matmul(_matmul(Fj, P), _transpose(Fj)), Q)
        # ---- update ----
        Hj = reverse_jacobian(h, xp)                    # dh/dx at predicted state
        hx = [v.value for v in h([Var(xi) for xi in xp])]
        y = [z[i] - hx[i] for i in range(len(z))]       # innovation
        S = _add(_matmul(_matmul(Hj, Pp), _transpose(Hj)), R)
        K = _matmul(_matmul(Pp, _transpose(Hj)), _inv(S))
        x = [xp[i] + v for i, v in enumerate(_matvec(K, y))]
        ImKH = _add(I, _matmul(K, Hj), sign=-1.0)
        # Joseph form keeps P symmetric positive-definite
        P = _add(_matmul(_matmul(ImKH, Pp), _transpose(ImKH)),
                 _matmul(_matmul(K, R), _transpose(K)))
        means.append(x[:])
        covs.append([r[:] for r in P])

    return {"filtered_means": means, "filtered_covariances": covs}
