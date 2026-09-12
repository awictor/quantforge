"""Optimal trade execution: the Almgren-Chriss model.

Liquidating ``X`` shares over ``T`` split into ``N`` intervals of length
``tau = T/N`` trades off market-impact cost against timing (price) risk. With a
linear temporary impact coefficient ``eta`` (per unit trading rate) and permanent
impact ``gamma``, and price volatility ``sigma``, the holdings trajectory that
minimizes ``E[cost] + lambda * Var[cost]`` is

    x_k = X * sinh(kappa (T - t_k)) / sinh(kappa T),   kappa^2 ~ lambda sigma^2 / eta~

so a risk-averse trader (``lambda > 0``) front-loads selling, while ``lambda = 0``
gives the linear (TWAP) schedule. This module computes the trajectory, the trade
list, the expected implementation-shortfall cost, its variance, and points on the
cost/variance efficient frontier. Pure standard library.
"""

import math


def _kappa(lam, sigma, eta, tau):
    """Urgency parameter kappa from the AC recursion.

    ``eta_tilde = eta - 0.5 gamma tau`` is folded into ``eta`` by the caller; here
    ``kappa_tilde^2 = lambda sigma^2 / eta``, ``cosh(kappa tau) = 1 + 0.5
    kappa_tilde^2 tau^2``, and kappa solves that. Returns 0 when ``lambda = 0``.
    """
    if lam <= 0.0:
        return 0.0
    kt2 = lam * sigma * sigma / eta
    ch = 1.0 + 0.5 * kt2 * tau * tau
    return math.acosh(ch) / tau


def execution_trajectory(total_shares, n_intervals, horizon, lam, sigma, eta):
    """Optimal holdings trajectory ``[x_0, x_1, ..., x_N]`` (Almgren-Chriss).

    ``x_k`` is the shares still held after interval ``k``. Starts at
    ``total_shares`` and ends at zero. For ``lam = 0`` the schedule is linear
    (TWAP); for ``lam > 0`` it is the ``sinh`` profile that liquidates faster
    early. ``eta`` is the temporary-impact coefficient.
    """
    if total_shares <= 0:
        raise ValueError("total_shares must be positive")
    if n_intervals < 1:
        raise ValueError("n_intervals must be a positive integer")
    if horizon <= 0 or eta <= 0 or sigma < 0:
        raise ValueError("horizon and eta must be positive, sigma non-negative")
    tau = horizon / n_intervals
    kappa = _kappa(lam, sigma, eta, tau)
    xs = []
    if kappa == 0.0:
        for k in range(n_intervals + 1):
            xs.append(total_shares * (1.0 - k / n_intervals))
        return xs
    sinh_kT = math.sinh(kappa * horizon)
    for k in range(n_intervals + 1):
        t_k = k * tau
        xs.append(total_shares * math.sinh(kappa * (horizon - t_k)) / sinh_kT)
    return xs


def execution_trades(trajectory):
    """Per-interval trade sizes ``n_k = x_{k-1} - x_k`` from a holdings trajectory.

    Positive sells that sum to the initial holdings.
    """
    return [trajectory[k - 1] - trajectory[k] for k in range(1, len(trajectory))]


def expected_cost(trajectory, horizon, gamma, eta):
    """Expected implementation-shortfall cost of a trajectory.

    Permanent impact contributes ``0.5 gamma X^2`` (independent of the path); the
    temporary impact contributes ``eta / tau * sum_k n_k^2`` for trades ``n_k``
    over intervals of length ``tau``. Returned in cash (price * shares) units.
    """
    n = len(trajectory) - 1
    if n < 1:
        raise ValueError("trajectory must have at least two points")
    tau = horizon / n
    trades = execution_trades(trajectory)
    X = trajectory[0]
    permanent = 0.5 * gamma * X * X
    temporary = (eta / tau) * sum(nk * nk for nk in trades)
    return permanent + temporary


def cost_variance(trajectory, horizon, sigma):
    """Timing-risk variance ``sigma^2 tau * sum_k x_k^2`` of the holdings path.

    The variance of execution cost from price moves while shares are still held;
    ``x_k`` are the holdings during each interval (using the end-of-interval
    holdings ``x_1..x_N``). Falls as liquidation is front-loaded.
    """
    n = len(trajectory) - 1
    if n < 1:
        raise ValueError("trajectory must have at least two points")
    tau = horizon / n
    return sigma * sigma * tau * sum(trajectory[k] ** 2 for k in range(1, n + 1))


def efficient_frontier_point(total_shares, n_intervals, horizon, lam, sigma,
                             eta, gamma):
    """One ``(expected_cost, variance)`` point for a given risk aversion ``lam``.

    Builds the :func:`execution_trajectory` at ``lam`` and returns its
    :func:`expected_cost` and :func:`cost_variance` -- sweeping ``lam`` traces the
    Almgren-Chriss efficient frontier (cost rises as variance falls).
    """
    traj = execution_trajectory(total_shares, n_intervals, horizon, lam, sigma, eta)
    return expected_cost(traj, horizon, gamma, eta), cost_variance(traj, horizon, sigma)
