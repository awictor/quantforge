"""Time-varying regression slope via a Kalman filter (dynamic hedge ratio).

A hedge ratio, factor loading, or pairs-trading beta is rarely constant. Model it
as a random walk observed through a regression:

    state:        beta_t = beta_{t-1} + w_t,     w_t ~ N(0, Q)
    observation:  y_t    = beta_t * x_t + v_t,   v_t ~ N(0, R)

This is a scalar Kalman filter with a *time-varying* observation loading ``x_t``:

    predict:  P_pred = P_{t-1} + Q
    gain:     K_t    = P_pred * x_t / (x_t^2 * P_pred + R)
    update:   beta_t = beta_pred + K_t (y_t - x_t * beta_pred)
              P_t    = (1 - K_t * x_t) * P_pred

With ``Q = 0`` and a diffuse prior the recursion is exactly recursive least
squares, so the final estimate equals the ordinary-least-squares slope
``sum(x y) / sum(x^2)``. A positive ``Q`` lets the slope drift, tracking a
changing relationship. Pure standard library.
"""


def kalman_regression_beta(x, y, process_var, obs_var, beta0=0.0, p0=1e6):
    """Filter a time-varying regression slope ``beta_t`` of ``y`` on ``x``.

    Parameters
    ----------
    x, y : sequence of float
        Regressor and response of equal length.
    process_var : float
        Slope random-walk variance ``Q`` (>= 0). ``0`` gives a static slope
        (recursive least squares); larger values let it drift faster.
    obs_var : float
        Observation noise variance ``R`` (> 0).
    beta0 : float
        Prior mean of the slope.
    p0 : float
        Prior variance of the slope. A large default gives a near-diffuse start so
        early data dominates.

    Returns
    -------
    (betas, variances) : (list[float], list[float])
        The filtered slope and its posterior variance at each step.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if process_var < 0:
        raise ValueError("process_var must be non-negative")
    if obs_var <= 0:
        raise ValueError("obs_var must be positive")

    n = len(x)
    beta = float(beta0)
    p = float(p0)
    betas = []
    variances = []
    for t in range(n):
        xt = float(x[t])
        # Predict.
        p_pred = p + process_var
        # Update with a time-varying loading x_t.
        denom = xt * xt * p_pred + obs_var
        k = p_pred * xt / denom
        beta = beta + k * (float(y[t]) - xt * beta)
        p = (1.0 - k * xt) * p_pred
        betas.append(beta)
        variances.append(p)
    return betas, variances
