"""Scalar Kalman filter for the local-level model.

The local-level (random-walk-plus-noise) state-space model is

    state:        x_t = x_{t-1} + w_t,      w_t ~ N(0, Q)
    observation:  y_t = x_t + v_t,          v_t ~ N(0, R)

a workhorse for tracking a slowly drifting level -- a time-varying mean, a dynamic
hedge ratio, or a smoothed signal. The Kalman filter is the exact
minimum-mean-square recursion for the hidden level ``x_t`` given observations up to
``t``:

    predict:  x_pred = x_{t-1},                 P_pred = P_{t-1} + Q
    gain:     K_t    = P_pred / (P_pred + R)
    update:   x_t    = x_pred + K_t (y_t - x_pred),   P_t = (1 - K_t) P_pred

The steady-state gain ``K`` solves the algebraic Riccati fixed point and depends
only on the signal-to-noise ratio ``q = Q / R``:

    K = (q + sqrt(q^2 + 4 q)) / (q + 2 + sqrt(q^2 + 4 q)).

As ``R -> 0`` the gain approaches 1 (trust each observation fully -> the filter
reproduces the data); as ``Q -> 0`` the gain decays like ``1/t`` and the estimate
converges to the running mean (recursive least squares). Pure standard library.
"""

import math


def kalman_local_level(observations, process_var, obs_var, x0=None, p0=None):
    """Run a scalar local-level Kalman filter over ``observations``.

    Parameters
    ----------
    observations : sequence of float
        The measured series ``y_t``.
    process_var : float
        State (process) noise variance ``Q`` (>= 0).
    obs_var : float
        Observation noise variance ``R`` (> 0).
    x0 : float, optional
        Prior mean for the level. Defaults to the first observation.
    p0 : float, optional
        Prior variance for the level. Defaults to ``obs_var`` (a diffuse-ish
        start). Larger values weight the data more at the outset.

    Returns
    -------
    (levels, variances, gains) : (list[float], list[float], list[float])
        The filtered level estimate, its posterior variance, and the Kalman gain
        at each step.
    """
    y = [float(v) for v in observations]
    n = len(y)
    if n == 0:
        return [], [], []
    if process_var < 0:
        raise ValueError("process_var must be non-negative")
    if obs_var <= 0:
        raise ValueError("obs_var must be positive")

    x = y[0] if x0 is None else float(x0)
    p = obs_var if p0 is None else float(p0)

    levels = []
    variances = []
    gains = []
    for t in range(n):
        # Predict.
        p_pred = p + process_var
        # Update.
        k = p_pred / (p_pred + obs_var)
        x = x + k * (y[t] - x)
        p = (1.0 - k) * p_pred
        levels.append(x)
        variances.append(p)
        gains.append(k)
    return levels, variances, gains


def kalman_steady_state_gain(process_var, obs_var):
    """Steady-state Kalman gain for the local-level model.

    Solves the algebraic Riccati fixed point in closed form; depends only on the
    signal-to-noise ratio ``q = process_var / obs_var``. Returns a gain in
    ``(0, 1)`` that rises toward 1 as the process noise dominates and toward 0 as
    the observation noise dominates.
    """
    if process_var < 0:
        raise ValueError("process_var must be non-negative")
    if obs_var <= 0:
        raise ValueError("obs_var must be positive")
    q = process_var / obs_var
    root = math.sqrt(q * q + 4.0 * q)
    return (q + root) / (q + 2.0 + root)
