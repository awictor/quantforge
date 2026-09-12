"""Ornstein-Uhlenbeck parameter estimation from a discrete series.

The OU process is the continuous-time mean-reverting SDE

    dX_t = kappa (theta - X_t) dt + sigma dW_t,

with mean-reversion speed ``kappa``, long-run level ``theta``, and diffusion
``sigma``. Sampled at a fixed step ``dt`` it is exactly a Gaussian AR(1):

    X_{t+1} = a + b X_t + eps,   b = exp(-kappa dt),
    a = theta (1 - b),           Var(eps) = sigma^2 (1 - b^2) / (2 kappa).

Fitting the AR(1) by least squares (its conditional MLE) and inverting these
relations recovers ``(kappa, theta, sigma)`` in closed form -- the standard
calibration behind a pairs / mean-reversion trade. The half-life of a deviation
is ``ln(2) / kappa``. Pure standard library.
"""

import math


def fit_ornstein_uhlenbeck(x, dt=1.0):
    """Estimate OU parameters ``(kappa, theta, sigma)`` from a sampled path.

    Parameters
    ----------
    x : sequence of float
        Observations sampled at uniform spacing ``dt``.
    dt : float
        Time between observations (in the same units as ``kappa`` is desired).

    Returns
    -------
    dict
        ``{"kappa", "theta", "sigma", "half_life"}``. ``kappa`` is the
        mean-reversion speed, ``theta`` the long-run mean, ``sigma`` the
        instantaneous volatility, and ``half_life = ln(2)/kappa``. Raises if the
        series is not mean-reverting (fitted ``b`` outside ``(0, 1)``).
    """
    n = len(x)
    if n < 3:
        raise ValueError("need at least 3 observations")
    if dt <= 0:
        raise ValueError("dt must be positive")

    xp = x[:-1]           # X_t
    xn = x[1:]            # X_{t+1}
    m = n - 1
    mx = sum(xp) / m
    my = sum(xn) / m
    sxx = sum((xp[i] - mx) ** 2 for i in range(m))
    if sxx <= 0.0:
        raise ValueError("regressor has zero variance")
    sxy = sum((xp[i] - mx) * (xn[i] - my) for i in range(m))

    b = sxy / sxx
    a = my - b * mx
    if not (0.0 < b < 1.0):
        raise ValueError("series is not mean-reverting (b outside (0, 1))")

    kappa = -math.log(b) / dt
    theta = a / (1.0 - b)

    # Residual variance of the AR(1) -> sigma via the exact OU relation.
    resid = [xn[i] - (a + b * xp[i]) for i in range(m)]
    resid_var = sum(r * r for r in resid) / m
    # Var(eps) = sigma^2 (1 - b^2) / (2 kappa)  ->  invert.
    sigma = math.sqrt(resid_var * 2.0 * kappa / (1.0 - b * b))

    return {
        "kappa": kappa,
        "theta": theta,
        "sigma": sigma,
        "half_life": math.log(2.0) / kappa,
    }
