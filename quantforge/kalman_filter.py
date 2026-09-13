"""General linear-Gaussian Kalman filter and RTS smoother.

State-space model
    x_t = F x_{t-1} + w_t,   w_t ~ N(0, Q)   (state transition)
    y_t = H x_t + v_t,        v_t ~ N(0, R)   (observation).

``kalman_filter`` runs the predict/update recursion, returning the filtered state
means and covariances and the data log-likelihood (from the one-step-ahead
prediction errors). ``kalman_smoother`` adds the Rauch-Tung-Striebel backward pass,
which uses future observations to refine each state, so its covariances never exceed
the filter's. Pure standard library on top of the matrix helpers.
"""

import math


def _matvec(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)]
            for i in range(n)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def _add(A, B, sign=1.0):
    return [[A[i][j] + sign * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _invert(A):
    from .portopt import _invert as inv
    return inv(A)


def kalman_filter(observations, F, H, Q, R, x0, P0):
    """Kalman filter over a sequence of observation vectors.

    ``F`` (state x state), ``H`` (obs x state), ``Q`` (state x state process cov),
    ``R`` (obs x obs measurement cov), ``x0``/``P0`` the initial state mean and
    covariance. Each entry of ``observations`` is an observation vector. Returns a
    dict with ``states`` (filtered means), ``covariances``, and ``log_likelihood``
    (the Gaussian log-likelihood of the observations from the prediction errors).
    """
    n = len(x0)
    x = list(x0)
    P = [row[:] for row in P0]
    states = []
    covs = []
    loglik = 0.0
    Ft = _transpose(F)
    Ht = _transpose(H)
    for y in observations:
        # Predict.
        xp = _matvec(F, x)
        Pp = _add(_matmul(_matmul(F, P), Ft), Q)
        # Innovation.
        yp = _matvec(H, xp)
        innov = [y[i] - yp[i] for i in range(len(y))]
        S = _add(_matmul(_matmul(H, Pp), Ht), R)
        Sinv = _invert(S)
        # Log-likelihood contribution.
        m = len(y)
        det = _det(S)
        quad = sum(innov[i] * sum(Sinv[i][j] * innov[j] for j in range(m))
                   for i in range(m))
        loglik += -0.5 * (m * math.log(2.0 * math.pi) + math.log(det) + quad)
        # Update (Kalman gain K = Pp H' Sinv).
        K = _matmul(_matmul(Pp, Ht), Sinv)
        x = [xp[i] + sum(K[i][j] * innov[j] for j in range(m)) for i in range(n)]
        KH = _matmul(K, H)
        ImKH = _add(_identity(n), KH, sign=-1.0)
        P = _matmul(ImKH, Pp)
        states.append(x[:])
        covs.append([row[:] for row in P])
    return {"states": states, "covariances": covs, "log_likelihood": loglik}


def _identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def _det(A):
    from .lu import determinant
    return determinant(A)


def kalman_smoother(observations, F, H, Q, R, x0, P0):
    """Rauch-Tung-Striebel smoother: filter then a backward refinement pass.

    Returns a dict with ``states`` (smoothed means) and ``covariances``. Each smoothed
    covariance is no larger (in the positive-definite sense; here checked on the
    trace) than the corresponding filtered covariance, since the smoother conditions
    on the whole series rather than only the past.
    """
    filt = kalman_filter(observations, F, H, Q, R, x0, P0)
    xf = filt["states"]
    Pf = filt["covariances"]
    T = len(xf)
    n = len(x0)
    Ft = _transpose(F)

    xs = [None] * T
    Ps = [None] * T
    xs[T - 1] = xf[T - 1][:]
    Ps[T - 1] = [row[:] for row in Pf[T - 1]]
    for t in range(T - 2, -1, -1):
        # Predicted covariance at t+1 from filtered state at t.
        Pp = _add(_matmul(_matmul(F, Pf[t]), Ft), Q)
        Ppinv = _invert(Pp)
        # Smoother gain C = Pf[t] F' Pp^{-1}.
        C = _matmul(_matmul(Pf[t], Ft), Ppinv)
        xpred = _matvec(F, xf[t])
        diff = [xs[t + 1][i] - xpred[i] for i in range(n)]
        xs[t] = [xf[t][i] + sum(C[i][j] * diff[j] for j in range(n)) for i in range(n)]
        Ps[t] = _add(Pf[t], _matmul(_matmul(C, _add(Ps[t + 1], Pp, sign=-1.0)),
                                    _transpose(C)))
    return {"states": xs, "covariances": Ps}
